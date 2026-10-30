"""ランの出力を 1 ステップ 1 個の ``.zarr.zip`` にまとめ、元と照合する (``r2d2plus-pack``)。

まとめるもの (``kinds``):

* ``qf``: 3D (``remap/qq``) → ``remap/qq/zarr/qq.NNNNNNNN.zarr.zip``
* ``qs``: slice → ``slice/zarr/<x|y|z>/qq.NNNNNNNN.zarr.zip``
* ``qt``: tau → ``tau/zarr/qq.NNNNNNNN.zarr.zip``

入れる変数は、元のファイルに書かれているものだけ (新出力形式では ``format.toml`` の
``variables`` などで、計算で出す ``pr``, ``te``, ``op`` は入れない)。従来形式は
ファイルにある変数を全部入れる。

照合 (元のファイルは 1 回だけ読む):

1. 元のファイルを読んだ配列を zarr ディレクトリに書き、読み戻して **ビット単位で** 比べる
2. zip にした後、全メンバーの CRC (``testzip``)、メンバーの名前と大きさが zarr
   ディレクトリと同じこと、配列の名前・shape・dtype・属性が同じことを見る

すでに zip があり元のファイルも残っているステップは、元と zip を読んでビット単位で
比べる (元を消すときの照合)。作業中の zarr と zip は ``.packing`` の付いた名前で作り、
照合が通ってから本来の名前にするので、読み手が途中のものを読むことはない。

元のファイルを消すのは ``delete_original=True`` のときだけで、さらに

* この呼び出しで照合が通ったステップ
* ランが走っていない (``data/param/run_summary.toml`` の status が ``normal_end`` か ``wall_limit``)
* ``wall_limit`` で止まったランでは、その種類の最後のステップではない (延長の前提なので残す)

のすべてを満たすものに限る。結果 (消したものを含む) は ``data/param/pack_report.toml`` に足していく。

``out`` を与えると、zip と記録はランではなく ``out/data/`` の下に作る (ランには何も書かない)。
このときは元のファイルを消さない。試しや、ランの外にまとめるとき用。
"""

import datetime
import json
import os
import re
import shutil
import socket
import time
import zipfile
from pathlib import Path

import numpy as np

try:  # Python 3.11+
    import tomllib
except ImportError:  # pragma: no cover
    tomllib = None

import pyR2D2
from pyR2D2.data_io import zarr_util

KINDS = ("qf", "qs", "qt")
REPORT = Path("param") / "pack_report.toml"
STOPPED = ("normal_end", "wall_limit")
TAUS = ["", "01", "001"]

CHECK = (
    "raw == zarr directory bit for bit (every array); "
    "zip: CRC of every member (testzip), member names and sizes == zarr directory, "
    "arrays (name, shape, dtype) and attrs"
)
CHECK_EXISTING = "raw == existing zip bit for bit (every array); zip: CRC of every member (testzip)"


def _same_bits(a, b):
    """単精度の配列がビット単位で等しいか (NaN の中身や -0.0 も区別する)。"""
    a = np.ascontiguousarray(a)
    b = np.ascontiguousarray(b)
    if a.shape != b.shape or a.dtype.kind != "f" or a.dtype.itemsize != 4:
        return False
    if b.dtype.kind != "f" or b.dtype.itemsize != 4:
        return False
    return np.array_equal(a.astype("=f4", copy=False).view(np.uint32), b.astype("=f4", copy=False).view(np.uint32))


def _toml_value(v):
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, (int, np.integer)):
        return str(int(v))
    if isinstance(v, float):
        return repr(v)
    if isinstance(v, (list, tuple)):
        return "[" + ", ".join(_toml_value(x) for x in v) + "]"
    return json.dumps(str(v), ensure_ascii=False)


def run_status(datadir):
    """``run_summary.toml`` の status (無ければ None)。"""
    path = Path(datadir) / "param" / "run_summary.toml"
    if not path.is_file() or tomllib is None:
        return None
    with open(path, "rb") as f:
        return tomllib.load(f).get("status", {}).get("status")


def parse_steps(text):
    """``"81,100-200"`` を数の集合にする。None は全部。"""
    if text is None:
        return None
    steps = set()
    for part in text.split(","):
        part = part.strip()
        if not part:
            continue
        if "-" in part:
            a, b = part.split("-", 1)
            steps.update(range(int(a), int(b) + 1))
        else:
            steps.add(int(part))
    return steps


class Unit:
    """まとめる単位 (種類・ステップ・slice の方向) と、その元のファイル・行き先。"""

    def __init__(self, kind, n, direc, raw, final):
        self.kind, self.n, self.direc = kind, n, direc
        self.raw = raw  # 元のファイル (読み手が読むもの全部)
        self.final = Path(final)  # 本来の zarr ディレクトリの名前 (zip はこれ + ".zip")

    @property
    def zip(self):
        return Path(str(self.final) + ".zip")

    @property
    def work_dir(self):
        return Path(str(self.final) + ".packing")

    @property
    def work_zip(self):
        return Path(str(self.final) + ".packing.zip")

    def label(self):
        return f"{self.kind}{'' if self.direc is None else ' ' + self.direc} n={self.n}"


class Packer:
    """
    1 つのランをまとめる

    Parameters
    ----------
    datadir : str
        ランの ``data/`` (またはランのディレクトリ)
    workers : int
        3D の元のファイルを読むスレッドの数
    log : callable
        進み具合を書く関数
    """

    def __init__(self, datadir, workers=1, log=print, out=None):
        datadir = Path(datadir)
        if (datadir / "data" / "param").is_dir():
            datadir = datadir / "data"
        self.datadir = datadir.resolve()
        self.outdir = None if out is None else (Path(out) / "data").resolve()
        run = self.datadir.parent
        if self.outdir is not None and (Path(out).resolve() == run or run in self.outdir.parents):
            raise ValueError("out must be outside the run")
        self.d = pyR2D2.Data(str(self.datadir), verbose=False)
        self.workers = workers
        self._log = log
        self._t0 = time.time()

    def dest(self, path):
        """ランの ``data/`` の中のパスを、書く場所 (``out`` があればその下) に読み替える。"""
        if self.outdir is None:
            return Path(path)
        return self.outdir / Path(path).relative_to(self.datadir)

    def log(self, *msg):
        self._log(f"[{time.time() - self._t0:8.1f}s]", *msg)

    def _compressed(self):
        return getattr(self.d.p, "output_format", "legacy") == "compressed"

    # ---- 何をまとめるか ----

    def steps(self, kind, direc=None):
        """元のファイルか zip がある ``kind`` のステップ (昇順)。"""
        dd = self.datadir
        found = set()
        if kind == "qf":
            for f in list((dd / "remap" / "qq").glob("*/*/qq.*")) + list((dd / "remap" / "qq").glob("qq.*")):
                parts = f.name.split(".")
                if len(parts) == 4 and parts[2].isdigit():
                    found.add(int(parts[2]))
            zdir = dd / "remap" / "qq" / "zarr"
        elif kind == "qs":
            for f in (dd / "slice").glob(f"qq{direc}*"):
                parts = f.name.split(".")
                if len(parts) == 4 and re.fullmatch(rf"qq{direc}(_ya?[ia]n)?", parts[0]) and parts[2].isdigit():
                    found.add(int(parts[2]))
            zdir = dd / "slice" / "zarr" / direc
        elif kind == "qt":
            for f in (dd / "tau").glob("qq.*"):
                parts = f.name.split(".")
                if len(parts) == 3 and parts[2].isdigit():
                    found.add(int(parts[2]))
            zdir = dd / "tau" / "zarr"
        else:
            raise ValueError(f"unknown kind {kind}")
        for f in self.dest(zdir).glob("qq.*.zarr.zip"):
            found.add(int(f.name.split(".")[1]))
        return sorted(found)

    def unit(self, kind, n, direc=None):
        d = self.d
        if kind == "qf":
            q = d.qf
            nps = []
            for ir0 in range(1, d.ixr + 1):
                for jr0 in range(1, d.jxr + 1):
                    np0 = d.np_ijr[ir0 - 1, jr0 - 1]
                    if ir0 == d.ir[np0] and jr0 == d.jr[np0]:
                        nps.append(np0)
            raw = [q._get_filepath_remap_qq(n, np0) for np0 in nps]
            return Unit(kind, n, None, raw, self.dest(q._get_filepath_remap_zarr(n)))
        if kind == "qs":
            q = d.qs
            raw = [
                q._get_filepath_slice(i, direc, n, pf)
                for i in range(len(q.xyz_slice_select(direc)))
                for pf in q._get_postfixes()
            ]
            return Unit(kind, n, direc, raw, self.dest(q._get_filepath_slice_zarr(n, direc)))
        q = d.qt
        return Unit(
            kind, n, None, [q._get_filepath_optical_depth(n)], self.dest(q._get_filepath_optical_depth_zarr(n))
        )

    def units(self, kinds=KINDS, steps=None):
        """まとめる単位と、種類 (と方向) ごとの最後のステップ。"""
        out = []
        for kind in kinds:
            direcs = [None]
            if kind == "qs":
                if not (self.datadir / "slice" / "params.dac").exists():
                    continue
                direcs = [c for c in "xyz" if len(self.d.qs.xyz_slice_select(c)) > 0]
            for direc in direcs:
                all_steps = self.steps(kind, direc)
                if not all_steps:
                    continue
                last = all_steps[-1]
                for n in all_steps:
                    if steps is None or n in steps:
                        out.append((self.unit(kind, n, direc), n == last))
        return out

    # ---- 元のファイルを読む ----

    def _keys_qf(self):
        if self._compressed():
            return list(self.d.p.format_info["variables"])
        return self.d.remap_keys + self.d.remap_keys_add

    def _read_qf(self, n, key):
        q = self.d.qf
        q.read(n=n, keys=key, max_workers=self.workers, zarr_flag=False)
        arr = q.__dict__[key]
        q.__dict__[key] = None
        return arr

    def _read_qs(self, unit):
        q = self.d.qs
        if self._compressed():
            keys = list(self.d.p.format_info["slice_variables"])
        else:
            q.read(n_slice=0, direc=unit.direc, n=unit.n)
            pf = q._get_postfixes()[0]
            keys = [
                k for k in self.d.remap_keys + self.d.remap_keys_add[:-1]
                if isinstance(q.__dict__.get(k + pf), np.ndarray)
            ]
        vars_dict, params, chunks = q.assemble(unit.n, unit.direc, keys)
        return vars_dict, params, chunks

    def _read_qt(self, unit):
        q = self.d.qt
        q.read(n=unit.n, zarr_flag=False)
        bases = list(q.value_keys) + ["height"] + [f"rt{g}" for g in range(1, 100)]
        vars_dict = {}
        for base in bases:
            # 新出力形式の pr, te は ro, se から計算したもの、he は height から作った倍精度
            if self._compressed() and base in ("pr", "te", "he"):
                continue
            for tau in TAUS:
                v = q.__dict__.get(base + tau)
                if isinstance(v, np.ndarray) and v.dtype.kind == "f" and v.dtype.itemsize == 4:
                    vars_dict[base + tau] = v
        chunk_max = 4096
        chunks = (1, min(self.d.jx, chunk_max), min(self.d.kx, chunk_max))
        return vars_dict, {}, chunks

    # ---- まとめる・照合する ----

    def _zip_and_check(self, unit, expected):
        """作業用の zarr ディレクトリを zip にし、照合して本来の名前にする。"""
        work_dir, work_zip = unit.work_dir, unit.work_zip
        zarr_util.zip_zarr(work_dir, zip_path=work_zip, overwrite=True)
        with zipfile.ZipFile(work_zip) as zf:
            bad = zf.testzip()
            members = {i.filename: i.file_size for i in zf.infolist()}
        if bad is not None:
            raise RuntimeError(f"CRC error in {work_zip}: {bad}")
        files = {
            f.relative_to(work_dir).as_posix(): f.stat().st_size
            for f in work_dir.rglob("*")
            if f.is_file()
        }
        if members != files:
            raise RuntimeError(f"members of {work_zip} differ from {work_dir}")
        root_dir = zarr_util.open_zarr_group(work_dir)
        root_zip = zarr_util.open_zarr_group(work_zip)
        if root_dir.attrs.asdict() != root_zip.attrs.asdict():
            raise RuntimeError(f"attrs of {work_zip} differ from {work_dir}")
        got = {k: (tuple(root_zip[k].shape), str(root_zip[k].dtype)) for k in root_zip.array_keys()}
        if got != expected:
            raise RuntimeError(f"arrays of {work_zip} differ: {got} != {expected}")
        os.replace(work_zip, unit.zip)
        shutil.rmtree(work_dir)

    def _save_checked(self, unit, vars_dict, params, chunks, mode):
        """配列を作業用の zarr に書き、読み戻してビット単位で比べる。"""
        chunks = tuple(int(c) for c in chunks)
        # tau の配列は 2D (jx, kx)。チャンクは OpticalDepth.compress と同じ (1, jx, kx) の後ろ 2 つ
        chunks2d = chunks[1:] if unit.kind == "qt" else None
        zarr_util.save(unit.work_dir, vars_dict, params, chunks2d=chunks2d, chunks3d=chunks, mode=mode)
        root = zarr_util.open_zarr_group(unit.work_dir)
        for key, arr in vars_dict.items():
            if key in ("x", "y", "z"):
                continue
            if not _same_bits(arr, root[key][...]):
                raise RuntimeError(f"{unit.label()}: {key} in the zarr differs from the binary files")
        return {k: (tuple(root[k].shape), str(root[k].dtype)) for k in vars_dict}

    def pack(self, unit):
        """元のファイルを 1 回読んで zip にまとめ、照合する。入れた変数の名前を返す。"""
        for p in (unit.work_dir, unit.work_zip):
            if p.is_dir():
                shutil.rmtree(p)
            elif p.exists():
                p.unlink()
        unit.final.parent.mkdir(parents=True, exist_ok=True)
        expected = {}
        if unit.kind == "qf":
            d = self.d
            chunks = d.qf.default_chunks3d(d.ix, d.jx, d.kx)
            params = {"i_start": 0, "i_size": d.ix, "j_start": 0, "j_size": d.jx, "k_start": 0, "k_size": d.kx}
            coords = {"x": d.x, "y": d.y, "z": d.z}
            expected.update(self._save_checked(unit, coords, params, chunks, "w"))
            keys = self._keys_qf()
            for key in keys:
                arr = self._read_qf(unit.n, key)
                expected.update(self._save_checked(unit, {key: arr}, None, chunks, "a"))
                del arr
        else:
            vars_dict, params, chunks = self._read_qs(unit) if unit.kind == "qs" else self._read_qt(unit)
            keys = list(vars_dict)
            expected.update(self._save_checked(unit, vars_dict, params, chunks, "w"))
            del vars_dict
        self._zip_and_check(unit, expected)
        return keys

    def verify_existing(self, unit):
        """すでにある zip を、元のファイルとビット単位で比べる。入っている変数の名前を返す。"""
        with zipfile.ZipFile(unit.zip) as zf:
            bad = zf.testzip()
        if bad is not None:
            raise RuntimeError(f"CRC error in {unit.zip}: {bad}")
        root = zarr_util.open_zarr_group(unit.zip)
        names = [k for k in root.array_keys() if k not in ("x", "y", "z")]
        if unit.kind == "qf":
            missing = [k for k in self._keys_qf() if k not in names]
            if missing:
                raise RuntimeError(f"{unit.zip} lacks {missing}")
            for key in self._keys_qf():
                if not _same_bits(self._read_qf(unit.n, key), root[key][...]):
                    raise RuntimeError(f"{unit.label()}: {key} in {unit.zip} differs from the binary files")
            return self._keys_qf()
        vars_dict = (self._read_qs(unit) if unit.kind == "qs" else self._read_qt(unit))[0]
        missing = [k for k in vars_dict if k not in names]
        if missing:
            raise RuntimeError(f"{unit.zip} lacks {missing}")
        for key, arr in vars_dict.items():
            if not _same_bits(arr, root[key][...]):
                raise RuntimeError(f"{unit.label()}: {key} in {unit.zip} differs from the binary files")
        return list(vars_dict)

    # ---- 記録 ----

    def record(self, entry):
        path = self.dest(self.datadir / REPORT)
        path.parent.mkdir(parents=True, exist_ok=True)
        new = not path.exists()
        with open(path, "a") as f:
            if new:
                f.write(
                    "# r2d2plus-pack (pyR2D2.pack) の記録。1 回の照合ごとに [[pack]] を 1 つ足す。\n"
                    "# result: packed (まとめて照合した) / verified (前からある zip を元と照合した) / failed\n"
                )
            f.write("\n[[pack]]\n")
            for key, value in entry.items():
                if value is None:
                    continue
                f.write(f"{key} = {_toml_value(value)}\n")

    # ---- 全体 ----

    def run(self, kinds=KINDS, steps=None, delete_original=False, dry_run=False):
        """
        まとめて照合し、許されていれば元のファイルを消す

        Returns
        -------
        bool
            照合に失敗したものが無ければ True
        """
        if delete_original and self.outdir is not None:
            raise ValueError("delete_original cannot be used with out (the zips are not in the run)")
        status = run_status(self.datadir)
        may_delete = delete_original and status in STOPPED
        # 走り終えたランは最後のステップも消してよい。wall_limit は延長の前提なので最後を残す
        keep_last = status != "normal_end"
        if delete_original and not may_delete:
            self.log(
                f"run status is {status!r} (param/run_summary.toml); "
                f"original files are not deleted unless it is one of {STOPPED}"
            )

        ok = True
        for unit, is_last in self.units(kinds, steps):
            present = [p for p in unit.raw if p.exists()]
            has_zip = unit.zip.exists()
            if not present:
                if not has_zip:
                    self.log(f"{unit.label()}: no files")
                continue
            if len(present) != len(unit.raw):
                self.log(f"{unit.label()}: only {len(present)}/{len(unit.raw)} binary files; skipped")
                continue
            if unit.final.is_dir():
                self.log(f"{unit.label()}: {unit.final} (zarr directory) exists; skipped")
                continue
            delete = may_delete and not (keep_last and is_last)
            action = "verify" if has_zip else "pack"
            if dry_run:
                note = "delete originals" if delete else ("keep (last step)" if may_delete else "keep")
                self.log(f"{unit.label()}: {action}, {len(present)} files, {note}")
                continue

            t0 = time.time()
            raw_bytes = sum(p.stat().st_size for p in present)
            entry = {
                "kind": unit.kind,
                "n": unit.n,
                "direc": unit.direc,
                "zip": unit.zip.relative_to(self.outdir or self.datadir).as_posix(),
                "raw_files": len(present),
                "raw_bytes": raw_bytes,
            }
            try:
                keys = self.pack(unit) if action == "pack" else self.verify_existing(unit)
                entry.update(result="packed" if action == "pack" else "verified", keys=keys)
                entry["zip_bytes"] = unit.zip.stat().st_size
                entry["check"] = CHECK if action == "pack" else CHECK_EXISTING
            except Exception as exc:  # 記録して次へ。元のファイルは消さない
                ok = False
                entry.update(result="failed", message=f"{type(exc).__name__}: {exc}")
                self.log(f"{unit.label()}: FAILED {exc}")
                delete = False

            deleted = []
            if delete:
                for p in present:
                    p.unlink()
                    deleted.append(p)
            entry["deleted_files"] = len(deleted)
            entry["deleted_bytes"] = raw_bytes if deleted else 0
            if deleted:
                parents = sorted({p.parent.relative_to(self.datadir).as_posix() for p in deleted})
                entry["deleted_example"] = deleted[0].relative_to(self.datadir).as_posix()
                entry["deleted_dirs"] = len(parents)
            if entry["result"] != "failed" and not deleted:
                why = "not requested"
                if delete_original:
                    why = f"run status {status!r}" if status not in STOPPED else "last step of a wall_limit run"
                entry["kept_because"] = why
            entry["seconds"] = round(time.time() - t0, 1)
            entry["date"] = datetime.datetime.now().astimezone().isoformat(timespec="seconds")
            entry["host"] = socket.gethostname()
            entry["pyR2D2"] = getattr(pyR2D2, "__version__", "")
            self.record(entry)
            if entry["result"] != "failed":
                self.log(
                    f"{unit.label()}: {entry['result']} {entry['zip_bytes'] / 1e9:.3f} GB"
                    f" ({raw_bytes / 1e9:.3f} GB raw), deleted {len(deleted)} files"
                )
        return ok
