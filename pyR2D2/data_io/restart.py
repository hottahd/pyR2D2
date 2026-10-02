"""
R2D2plus の restart を読む (R2D2plus DEC-074〜078・430・557、specification 6.1・11.34)。

restart は通常出力 (legacy / compressed) とは別の専用形式で、倍精度・C++ の軸順 (**z が鉛直**)・
margin 無しの物理セルだけを持つ。ここでは書かない (書き換えは R2D2plus の ``r2d2plus_upgrade`` の役目)。

配置::

    data/restart/<slot>/meta.toml                         slot は "e" / "o"(交代)と "00000010"(恒久)
    data/restart/<slot>/state.<rank>.bin                  v1 (format_version = 1)
    data/restart/<slot>/<generation>/state.<rank>.bin     v2 (meta.toml の [restart] generation)
    data/restart/<slot>/<generation>/{yin,yang}/state.<rank>.bin   YinYang
    data/restart/<slot>/<generation>/extra.<name>.bin     全ランク共有の補助 (AFC など)
    data/restart/<slot>/<generation>/[<panel>/]extra.<name>.<rank>.bin   ランクごとの補助

各 state ファイルは ro, vx, vy, vz, bx, by, bz, se, ps の 9 変数を変数の順に、ランクの物理セルを
x → y → z (z が最速) で並べた little endian の float64。``[restart] compression_level`` が 0 でなければ
``r2d2plus-z`` の変数 "state" に同じ並びが入る。ランク番号は MPI の直交トポロジの順 (z が最速)。
meta.toml の ``[checksum] values`` は保存したバイト列 (圧縮後) の XXH3-64 で、world rank 順
(YinYang は Yin 0..P-1、Yang P..2P-1)。
"""

from dataclasses import dataclass, field
from pathlib import Path
import tomllib

import numpy as np

from .compressed import RESTART_COMPONENTS, _xxh3_hex, read_restart_state

__all__ = ["RESTART_VARIABLES", "Restart", "read_restart", "restart_slots"]

#: 9 変数の名前と順番 (:data:`pyR2D2.data_io.compressed.RESTART_COMPONENTS` と同じ)
RESTART_VARIABLES = tuple(RESTART_COMPONENTS)


@dataclass
class Restart:
    """:func:`read_restart` の結果。"""

    path: Path
    meta: dict
    #: ``{name: ndarray (nx, ny, nz)}``。YinYang では ``{"yin": {...}, "yang": {...}}``
    state: dict
    #: 補助 payload の名前 → ファイルの一覧 (読むのは :meth:`extra_bytes`)
    extras: dict = field(default_factory=dict)

    @property
    def time(self):
        return float(self.meta["restart"]["time"])

    @property
    def step(self):
        return int(self.meta["restart"]["step"])

    @property
    def yinyang(self):
        return bool(self.meta.get("mpi", {}).get("yinyang", False))

    @property
    def shape(self):
        g = self.meta["grid"]
        return (int(g["nx"]), int(g["ny"]), int(g["nz"]))

    def extra_bytes(self, name, rank=0):
        """補助 payload の生のバイト列。形式は payload ごとに R2D2plus が決める。"""
        files = self.extras[name]
        return Path(files[0] if len(files) == 1 else files[rank]).read_bytes()


def _restart_dir(path):
    path = Path(path)
    if (path / "meta.toml").exists():
        return path.parent, path.name
    if (path / "restart").is_dir():
        path = path / "restart"
    elif (path / "data" / "restart").is_dir():
        path = path / "data" / "restart"
    return path, None


def restart_slots(path):
    """
    restart の slot を ``{slot: (step, time)}`` で返す。

    Parameters
    ----------
    path : str or pathlib.Path
        ランのディレクトリ、``data/``、``data/restart/`` のどれか
    """
    root, _ = _restart_dir(path)
    out = {}
    for slot in sorted(p for p in root.iterdir() if (p / "meta.toml").exists()):
        meta = tomllib.loads((slot / "meta.toml").read_text())
        out[slot.name] = (int(meta["restart"]["step"]), float(meta["restart"]["time"]))
    return out


def read_restart(path, slot="latest", verify=None):
    """
    R2D2plus の restart を読み、全体の配列に組み立てる。

    Parameters
    ----------
    path : str or pathlib.Path
        ランのディレクトリ、``data/``、``data/restart/``、または slot のディレクトリ
        (``meta.toml`` があるところ)
    slot : str, optional
        ``"latest"`` (既定。step が最大の slot)、``"e"``、``"o"``、``"00000010"`` など。
        ``path`` が slot のディレクトリなら無視する。
    verify : bool or None, optional
        meta.toml の XXH3-64 と照合するか。``None`` (既定) なら ``xxhash`` モジュールが
        あるときだけ照合する。``True`` で ``xxhash`` が無ければ例外。

    Returns
    -------
    Restart
        ``state[name]`` は shape ``(nx, ny, nz)`` の float64 (C++ の軸順、**z が鉛直**。
        legacy 出力の Fortran 順 (鉛直が先) とは違う)。YinYang は ``state["yin"][name]`` など。

    Examples
    --------
    .. code-block:: python

        r = pyR2D2.read_restart("run/d001")          # 最新の slot
        r.time, r.step, r.shape
        r.state["ro"][:, :, -1]                      # 上端の密度
    """
    root, given = _restart_dir(path)
    if given is None:
        slots = restart_slots(root)
        if not slots:
            raise FileNotFoundError(f"{root}: no restart slot with meta.toml")
        if slot == "latest":
            given = max(slots, key=lambda s: (slots[s][0], slots[s][1]))
        elif slot not in slots:
            raise FileNotFoundError(f"{root}: no slot {slot!r} (have {sorted(slots)})")
        else:
            given = slot
    slot_dir = root / given
    meta = tomllib.loads((slot_dir / "meta.toml").read_text())
    rst = meta["restart"]
    version = rst.get("format_version")
    if version == 2:
        payload = slot_dir / rst["generation"]
    elif version == 1:
        payload = slot_dir
    else:
        raise ValueError(f"{slot_dir}: unsupported restart format_version {version!r}")

    if verify and _xxh3_hex(b"") is None:
        raise ImportError("verify=True needs the xxhash module (`pip install xxhash`)")

    g, m = meta["grid"], meta["mpi"]
    nx, ny, nz = int(g["nx"]), int(g["ny"]), int(g["nz"])
    sx, sy, sz = int(m["size_x"]), int(m["size_y"]), int(m["size_z"])
    if nx % sx or ny % sy or nz % sz:
        raise ValueError(f"{slot_dir}: grid {nx, ny, nz} is not divisible by mpi {sx, sy, sz}")
    lx, ly, lz = nx // sx, ny // sy, nz // sz

    panels = ["yin", "yang"] if m.get("yinyang") else [None]
    nranks = sx * sy * sz
    state = {}
    for panel in panels:
        full = {name: np.empty((nx, ny, nz)) for name in RESTART_VARIABLES}
        for rank in range(nranks):
            # 1 ランク分の読み込み・照合・復号は read_restart_state に任せる (実装を 1 つにする)
            local, _ = read_restart_state(slot_dir, rank, state_subdir=panel or "", verify=verify)
            ix, rem = divmod(rank, sy * sz)
            iy, iz = divmod(rem, sz)
            for name in RESTART_VARIABLES:
                full[name][ix * lx:(ix + 1) * lx, iy * ly:(iy + 1) * ly, iz * lz:(iz + 1) * lz] = local[name]
        state[panel] = full
    if panels == [None]:
        state = state[None]

    extras = {}
    for name, spec in meta.get("extras", {}).items():
        if spec.get("shared"):
            extras[name] = [payload / f"extra.{name}.bin"]
        else:
            extras[name] = [(payload / panel if panel else payload) / f"extra.{name}.{r:08d}.bin"
                            for panel in panels for r in range(nranks)]
    return Restart(path=slot_dir, meta=meta, state=state, extras=extras)
