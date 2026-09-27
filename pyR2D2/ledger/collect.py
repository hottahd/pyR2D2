"""ランディレクトリ ``<project>/run/dNNN`` から台帳の管理列の値を集める。

:py:class:`pyR2D2.Data` は使わない (remap 情報や EOS 表まで読むので重い)。
読むのは次の小さなファイルだけである。

* ``data/param/params.dac`` (Fortran 版と R2D2plus の従来形式が書く。
  C++ の 2D ランには無い)、``back.dac`` (座標と RSST の xi)
* ``data/param/run_summary.toml`` (R2D2plus DEC-578。状態・時刻・commit など)
* ``data/param/origin.toml`` (起点。Fortran 版は ``data/cont_log.txt``)
* ``data/param/runs/<最新>/effective_config.txt`` と ``run.log`` の末尾
* ``data/param/format.toml``、``data/restart/<slot>/meta.toml`` の ``[modifiers]``
* ``data/param/ledger_sent.toml`` (前回送ったときの要約のハッシュ)

容量 (``data/`` の使用量) だけは全ファイルの stat が要るので、``du=False`` で省ける。
"""

import datetime
import hashlib
import json
import os
import re
import socket
import tomllib
import zipfile
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from ..constant import constant
from .units import Quantity

# 管理列 (この順でシートに無い列を足す)。見出しは単位を除いた名前で探す。
# (名前, 別名, 量の列か)
COLUMNS = [
    ("Case ID", (), False),
    ("Mstar", (), False),
    ("(ix,jx,kx)", (), False),
    ("xmin", (), True),
    ("xmax", (), True),
    ("ymin", (), True),
    ("ymax", (), True),
    ("zmin", (), True),
    ("zmax", (), True),
    ("uni", (), False),
    ("dx", (), True),
    ("m ray", (), False),
    ("dtout", (), True),
    ("dtout_tau", (), True),
    ("al", (), False),
    ("RSST", (), False),
    ("Om", (), False),
    ("Geometry", ("Gemetry",), False),
    ("origin", (), False),
    ("update time", (), False),
    ("Server", (), False),
    ("状態", (), False),
    ("到達時刻", (), True),
    ("t_end", (), True),
    ("step", (), False),
    ("形式", (), False),
    ("容量", (), True),
    ("加工", (), False),
    ("commit", (), False),
    ("backend", (), False),
    ("ranks", (), False),
    ("ms/step", (), False),
]
COLUMN_NAMES = [c[0] for c in COLUMNS]
QUANTITY_COLUMNS = {c[0] for c in COLUMNS if c[2]}

# 単位を共有する量の列の組 (値を合わせて 1 つの単位を選び、変えるときは一緒に変える)。
# dtout と dtout_tau は独立。dx の「下端 上端」の 2 つの数はもとから同じ単位。
UNIT_GROUPS = [
    ("xmin", "xmax"),
    ("ymin", "ymax"),
    ("zmin", "zmax"),
    ("到達時刻", "t_end"),
]

# ハッシュに入れない列 (容量は du が要るので status では測らない)
HASH_EXCLUDE = {"容量"}

STATUS_LABELS = {
    "running": "走行中",
    "normal_end": "正常終了",
    "wall_limit": "時間切れ停止",
    "error": "異常終了",
}
STALE_LABEL = "不明"
STALE_AFTER = datetime.timedelta(days=1)

SENT_FILE = "ledger_sent.toml"

_CASE_RE = re.compile(r"^d(\d+)$")


def default_server():
    """既定のサーバ名 (ホスト名の最初の区切りまで)。"""
    return socket.gethostname().split(".")[0]


def case_row(caseid):
    """Case ID ``dNNN`` を置くシートの行番号 (1 始まり。d001 → 2)。dNNN でなければ ``None``。"""
    m = _CASE_RE.match(caseid)
    if not m:
        return None
    return int(m.group(1)) + 1


# ---------------------------------------------------------------------------
# 小さなファイルの読み込み
# ---------------------------------------------------------------------------


def _read_toml(path):
    try:
        with open(path, "rb") as f:
            return tomllib.load(f)
    except (OSError, tomllib.TOMLDecodeError):
        return None


def read_params_dac(path):
    """``params.dac`` ("値 名前 型" の行) を dict で返す。Parameters と同じ読み方。"""
    out = {}
    with open(path, "r", errors="replace") as f:
        for line in f:
            parts = line.split()
            if len(parts) == 2:
                parts = [""] + parts
            elif len(parts) < 2:
                continue
            value, key, vtype = parts[0], parts[1], parts[2]
            try:
                if vtype == "i":
                    out[key] = int(value)
                elif vtype == "d":
                    out[key] = float(value)
                elif vtype == "c":
                    out[key] = value
                elif vtype == "l":
                    out[key] = value != "F"
            except ValueError:
                continue
    return out


def read_back_x_xi(param_dir, p):
    """``back.dac`` から鉛直座標 x と xi (のりしろを除いた内部) を読む。

    読めなければ ``(None, None)``。
    """
    try:
        endian = "<" if p.get("swap", 0) == 0 else ">"
        margin = p["margin"]
        ixg = p["ix0"] * p["nx"] + 2 * margin * (p["xdcheck"] - 1)
        jxg = p["jx0"] * p["ny"] + 2 * margin * (p["ydcheck"] - 1)
        kxg = p["kx0"] * p["nz"] + 2 * margin * (p["zdcheck"] - 1)
        # Parameters と同じ並び: head, x, y, z, 1 次元量 20 本 (最後が xi), tail
        names = [
            "pr0", "te0", "ro0", "se0", "en0", "op0", "tu0", "dsedr0", "dtedr0",
            "dprdro", "dprdse", "dtedro", "dtedse", "dendro", "dendse", "gx",
            "cp", "fa", "sa", "xi",
        ]  # fmt: skip
        dtype = np.dtype(
            [("head", endian + "i"), ("x", endian + "d", (ixg,)),
             ("y", endian + "d", (jxg,)), ("z", endian + "d", (kxg,))]
            + [(n, endian + "d", (ixg,)) for n in names]
            + [("tail", endian + "i")]
        )  # fmt: skip
        back = np.fromfile(param_dir / "back.dac", dtype=dtype, count=1)
        if back.size != 1:
            return None, None
        x = back["x"][0][margin : ixg - margin]
        xi = back["xi"][0][margin : ixg - margin]
        return x, xi
    except (OSError, KeyError, ValueError, TypeError):
        return None, None


def read_effective_config(path):
    """``effective_config.txt`` ("key = value" の行) を dict で返す。数は float にする。"""
    out = {}
    try:
        with open(path, "r", errors="replace") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#") or "=" not in line:
                    continue
                key, _, value = line.partition("=")
                value = value.strip()
                if len(value) >= 2 and value[0] == value[-1] == '"':
                    value = value[1:-1]
                else:
                    try:
                        value = float(value)
                    except ValueError:
                        if value in ("true", "false"):
                            value = value == "true"
                out[key.strip()] = value
    except OSError:
        pass
    return out


def latest_segment(param_dir):
    """``param/runs/<区間>`` のうち名前 (起動時刻) が最新のもの。無ければ ``None``。"""
    runs = param_dir / "runs"
    try:
        segs = sorted(d for d in runs.iterdir() if d.is_dir())
    except OSError:
        return None
    return segs[-1] if segs else None


_MS_RE = re.compile(r"time loop .*?\(([0-9.eE+-]+) ms/step\)")


def read_ms_per_step(segment):
    """最新区間の ``run.log`` の末尾から ``time loop ... (X ms/step)`` を読む。"""
    if segment is None:
        return None
    path = segment / "run.log"
    try:
        with open(path, "rb") as f:
            f.seek(0, os.SEEK_END)
            size = f.tell()
            f.seek(max(0, size - 65536))
            tail = f.read().decode("utf-8", errors="replace")
    except OSError:
        return None
    found = _MS_RE.findall(tail)
    if not found:
        return None
    try:
        return float(found[-1])
    except ValueError:
        return None


def latest_restart_meta(data_dir):
    """``restart/<slot>/meta.toml`` のうち step が最大のものを返す。"""
    best, best_step = None, None
    try:
        slots = [d for d in (data_dir / "restart").iterdir() if d.is_dir()]
    except OSError:
        return None
    for slot in slots:
        meta = _read_toml(slot / "meta.toml")
        if not meta:
            continue
        step = meta.get("restart", {}).get("step", -1)
        if best_step is None or step > best_step:
            best, best_step = meta, step
    return best


def disk_usage(path):
    """``du`` 相当 (割り当てブロックの合計、バイト)。ハードリンクは 1 回だけ数える。"""
    total = 0
    seen = set()
    stack = [str(path)]
    while stack:
        top = stack.pop()
        try:
            with os.scandir(top) as it:
                for e in it:
                    try:
                        st = e.stat(follow_symlinks=False)
                    except OSError:
                        continue
                    if e.is_dir(follow_symlinks=False):
                        stack.append(e.path)
                    if st.st_nlink > 1:
                        key = (st.st_dev, st.st_ino)
                        if key in seen:
                            continue
                        seen.add(key)
                    total += getattr(st, "st_blocks", 0) * 512 or st.st_size
        except OSError:
            continue
    return total


def _fortran_time(data_dir, p):
    """Fortran 版のランの到達時刻 (``time/mhd/t.dac.<nd>``)。読めなければ ``None``。"""
    try:
        with open(data_dir / "param" / "nd.dac") as f:
            nd = int(f.read().split()[0])
        endian = "<" if p.get("swap", 0) == 0 else ">"
        name = f"t.dac.{nd:08d}"
        path = data_dir / "time" / "mhd" / name
        if path.exists():
            return float(np.fromfile(path, endian + "d", 1)[0])
        zpath = data_dir / "time.zip"
        if zpath.exists():
            with zipfile.ZipFile(zpath) as zf:
                buf = zf.read(f"time/mhd/{name}")
            return float(np.frombuffer(buf, endian + "d", 1)[0])
    except (OSError, ValueError, IndexError, KeyError):
        pass
    return None


def _parse_updated(text):
    try:
        return datetime.datetime.strptime(text, "%Y-%m-%dT%H:%M:%S%z")
    except (TypeError, ValueError):
        try:
            return datetime.datetime.fromisoformat(text)
        except (TypeError, ValueError):
            return None


# ---------------------------------------------------------------------------
# 1 本のラン
# ---------------------------------------------------------------------------


@dataclass
class RunRecord:
    """1 本のランの管理列の値。値は文字列、:class:`Quantity`、または ``None`` (空欄)。"""

    caseid: str
    run_dir: Path
    values: dict = field(default_factory=dict)
    run_summary_sha256: str = ""
    fortran: bool = False

    def digest(self):
        """管理列の値 (容量を除く) のハッシュ。送信済みかの判定に使う。"""
        canon = {}
        for name in COLUMN_NAMES:
            if name in HASH_EXCLUDE:
                continue
            v = self.values.get(name)
            canon[name] = v.canonical() if isinstance(v, Quantity) else v
        blob = json.dumps(canon, ensure_ascii=False, sort_keys=True)
        return hashlib.sha256(blob.encode("utf-8")).hexdigest()[:16]

    @property
    def sent_path(self):
        return self.run_dir / "data" / "param" / SENT_FILE

    def read_sent(self):
        return _read_toml(self.sent_path)

    def unsent(self):
        """``"never"`` (一度も送っていない)、``True`` (変更あり)、``False``。"""
        sent = self.read_sent()
        if not sent or "sent" not in sent:
            return "never"
        return sent["sent"].get("hash") != self.digest()


def _run_summary_sha(param_dir):
    try:
        return hashlib.sha256((param_dir / "run_summary.toml").read_bytes()).hexdigest()
    except OSError:
        return ""


def collect_run(run_dir, server=None, du=True, now=None):
    """ランディレクトリ 1 本から :class:`RunRecord` を作る。"""
    run_dir = Path(run_dir)
    data_dir = run_dir / "data"
    param_dir = data_dir / "param"
    server = server or default_server()
    now = now or datetime.datetime.now(datetime.timezone.utc)

    p = None
    if (param_dir / "params.dac").exists():
        p = read_params_dac(param_dir / "params.dac")
    rs = _read_toml(param_dir / "run_summary.toml")
    st = (rs or {}).get("status", {})
    rn = (rs or {}).get("run", {})
    origin = _read_toml(param_dir / "origin.toml")
    fmt = _read_toml(param_dir / "format.toml")
    segment = latest_segment(param_dir)
    eff = read_effective_config(segment / "effective_config.txt") if segment else {}

    rec = RunRecord(caseid=run_dir.name, run_dir=run_dir, fortran=rs is None)
    rec.run_summary_sha256 = _run_summary_sha(param_dir)
    v = rec.values
    v["Case ID"] = rec.caseid

    # 幾何
    geometry = None
    if p is not None and p.get("geometry"):
        geometry = p["geometry"]
    elif rn:
        if rn.get("yinyang"):
            geometry = "YinYang"
        else:
            g = str(rn.get("geometry", ""))
            geometry = {"cartesian": "Cartesian", "spherical": "Spherical"}.get(
                g.lower(), g or None
            )
    v["Geometry"] = geometry

    # 恒星と格子
    rstar = None
    if p is not None:
        rstar = p.get("rstar", constant.RSUN)
        mstar = p.get("mstar")
        v["Mstar"] = f"{mstar / constant.MSUN:.2f}" if mstar else "1.00"
        try:
            ix = p["nx"] * p["ix0"]
            jx = p["ny"] * p["jx0"]
            kx = p["nz"] * p["kx0"]
            if geometry == "YinYang":  # Parameters と同じ球座標への換算後の格子数
                jx = jx * 2
                kx = jx * 2
            v["(ix,jx,kx)"] = f"{ix} {jx} {kx}"
        except KeyError:
            ix = None
    elif rn and all(k in rn for k in ("nx", "ny", "nz")):
        # C++ の命名 (z が鉛直) を Fortran の (x 鉛直, y, z) に並べ替える
        v["(ix,jx,kx)"] = f"{rn['nz']} {rn['nx']} {rn['ny']}"
    if rstar is None and fmt and "rstar" in fmt:
        rstar = float(fmt["rstar"])

    # 領域
    if p is not None and "xmin" in p:
        dom = {k: p.get(k) for k in ("xmin", "xmax", "ymin", "ymax", "zmin", "zmax")}
    elif geometry == "Cartesian" and "domain.zmin" in eff:
        # params.dac の無い C++ ラン: C++ の z (鉛直) → x、x → y、y → z
        dom = {
            "xmin": eff.get("domain.zmin"), "xmax": eff.get("domain.zmax"),
            "ymin": eff.get("domain.xmin"), "ymax": eff.get("domain.xmax"),
            "zmin": eff.get("domain.ymin"), "zmax": eff.get("domain.ymax"),
        }  # fmt: skip
    else:
        dom = {}
    for key in ("xmin", "xmax"):
        x = dom.get(key)
        if not isinstance(x, float) or rstar is None:
            continue
        if geometry in ("Spherical", "YinYang"):
            v[key] = Quantity("ratio", (x / rstar,))
        else:
            v[key] = Quantity("length", (x - rstar,))
    if geometry == "YinYang":
        for key, deg in (("ymin", 0.0), ("ymax", 180.0), ("zmin", -180.0), ("zmax", 180.0)):
            v[key] = Quantity("angle", (deg,))
    else:
        for key in ("ymin", "ymax", "zmin", "zmax"):
            y = dom.get(key)
            if not isinstance(y, float):
                continue
            if geometry == "Spherical":
                v[key] = Quantity("angle", (float(np.degrees(y)),))
            else:
                v[key] = Quantity("length", (y,))

    # 格子幅・RSST (back.dac)
    if p is not None:
        v["uni"] = "F" if p.get("ununiform_flag") else "T"
        x, xi = read_back_x_xi(param_dir, p)
        if x is not None and x.size >= 2:
            v["dx"] = Quantity("length", (float(x[1] - x[0]), float(x[-1] - x[-2])))
        if xi is not None and xi.size:
            v["RSST"] = "F" if float(xi.max()) == 1.0 else "T"
        rte = p.get("rte")
        if isinstance(rte, bool):
            ngroups = p.get("rte_num_groups", 1)
            if not rte:
                v["m ray"] = "F"
            elif ngroups and ngroups > 1:
                v["m ray"] = f"T ({ngroups}群)"
            else:
                v["m ray"] = "T"
        elif rte is not None:
            v["m ray"] = str(rte)
        if "potential_alpha" in p:
            v["al"] = f"{p['potential_alpha']:.2f}"
        if "omfac" in p:
            v["Om"] = f"{p['omfac']:.1f}"

    # 出力の間隔と終了時刻
    def _time(value):
        if isinstance(value, (int, float)) and not isinstance(value, bool):
            return Quantity("time", (float(value),))
        return None

    dtout = p.get("dtout") if p is not None else eff.get("output.dtout")
    dtout_tau = p.get("dtout_tau") if p is not None else eff.get("output.dtout_tau")
    v["dtout"] = _time(dtout)
    v["dtout_tau"] = _time(dtout_tau)
    t_end = eff.get("time.t_end", rn.get("t_end") if rn else None)
    if t_end is None and p is not None:
        t_end = p.get("tend")
    v["t_end"] = _time(t_end) if t_end else None

    # 起点
    if origin and "origin" in origin:
        o = origin["origin"]
        name = Path(str(o.get("parent", ""))).name or "?"
        slot = o.get("slot")
        v["origin"] = f"{name} ({slot})" if slot else name
    elif (data_dir / "cont_log.txt").exists():
        try:
            with open(data_dir / "cont_log.txt") as f:
                v["origin"] = f.readlines()[6][-11:-7]
        except (OSError, IndexError):
            v["origin"] = "N/A"
    else:
        v["origin"] = "N/A"

    # 状態・時刻
    v["Server"] = server
    if rs is not None:
        updated = _parse_updated(st.get("updated"))
        if updated is not None:
            v["update time"] = updated.strftime("%Y-%m-%d %H:%M:%S")
        code = st.get("status")
        label = STATUS_LABELS.get(code, code)
        if code == "running" and updated is not None:
            if updated.tzinfo is None:
                updated = updated.astimezone()
            if now - updated > STALE_AFTER:
                label = STALE_LABEL
        v["状態"] = label
        if "time" in st:
            v["到達時刻"] = _time(float(st["time"]))
        if "step" in st:
            v["step"] = str(st["step"])
        commit = st.get("commit")
        if commit:
            v["commit"] = commit[:8] + ("-dirty" if st.get("dirty") else "")
        if st.get("backend"):
            v["backend"] = str(st["backend"])
        if "world_size" in st:
            v["ranks"] = str(st["world_size"])
    else:
        if p is not None:
            t = _fortran_time(data_dir, p)
            v["到達時刻"] = _time(t) if t is not None else None
            if "npe" in p:
                v["ranks"] = str(p["npe"])
        try:
            mtime = (param_dir / "nd.dac").stat().st_mtime
            v["update time"] = datetime.datetime.fromtimestamp(mtime).strftime(
                "%Y-%m-%d %H:%M:%S"
            )
        except OSError:
            pass

    # 出力形式
    fmt_name = eff.get("output.format")
    if not fmt_name and fmt:
        fmt_name = fmt.get("format")
    if not fmt_name and (p is not None or rs is not None):
        fmt_name = "legacy"
    v["形式"] = fmt_name or None

    # 施した加工 (最新のリスタートの [modifiers] applied)
    meta = latest_restart_meta(data_dir)
    if meta:
        applied = meta.get("modifiers", {}).get("applied", [])
        v["加工"] = "; ".join(str(a) for a in applied) if applied else None

    ms = read_ms_per_step(segment)
    v["ms/step"] = f"{ms:.2f}" if ms is not None else None

    if du and data_dir.is_dir():
        v["容量"] = Quantity("bytes", (float(disk_usage(data_dir)),))

    for name in COLUMN_NAMES:
        if name == "容量" and not du:
            continue  # 測っていない列は書かない (シートの値を空で上書きしない)
        v.setdefault(name, None)
    return rec


# ---------------------------------------------------------------------------
# プロジェクト
# ---------------------------------------------------------------------------


def _sort_key(path):
    m = _CASE_RE.match(path.name)
    return (0, int(m.group(1)), "") if m else (1, 0, path.name)


def find_runs(path):
    """``<project>/run`` の下のランディレクトリ (``data/`` を持つもの) を返す。

    ``path`` 自体が ``data/`` を持てば、そのラン 1 本とみなす。
    """
    path = Path(path)
    if (path / "data").is_dir():
        return [path]
    runs = [d for d in path.iterdir() if d.is_dir() and (d / "data").is_dir()]
    return sorted(runs, key=_sort_key)


def logical_path(path):
    """シンボリックリンクを**たどらない**絶対パス。

    ``run/d003/data -> /storage/d003_output`` のように data を別の場所へ置く運用では、
    ``resolve()`` するとラン(d003)ではなく置き場所を指してしまう(Codex の指摘、2026-09-27)。
    """
    return Path(os.path.abspath(os.fspath(path)))


def run_dir_of(datadir):
    """``<run>/data`` から ``<run>`` を返す(data がシンボリックリンクでもランを指す)。"""
    return logical_path(datadir).parent


def project_name(path):
    """シート名の既定 (プロジェクトのディレクトリ名)。"""
    path = logical_path(path)
    if (path / "data").is_dir():  # ラン 1 本: <project>/run/dNNN
        return path.parents[1].name
    return path.parent.name


def write_sent(record, sheet):
    """送信の記録 ``data/param/ledger_sent.toml`` を書く。"""
    now = datetime.datetime.now().astimezone().strftime("%Y-%m-%dT%H:%M:%S%z")

    def q(s):
        return json.dumps(str(s), ensure_ascii=False)

    text = (
        "# r2d2plus-ledger push で台帳へ送った要約 (管理列のハッシュ)。\n"
        "# 手で消すと、次の status で「未送信」と表示される。\n"
        "[sent]\n"
        f"hash = {q(record.digest())}\n"
        f"run_summary_sha256 = {q(record.run_summary_sha256)}\n"
        f"time = {q(now)}\n"
        f"sheet = {q(sheet)}\n"
        f"server = {q(record.values.get('Server') or '')}\n"
    )
    path = record.sent_path
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(text, encoding="utf-8")
    os.replace(tmp, path)
