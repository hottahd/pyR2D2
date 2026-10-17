"""
R2D2plus の光球スラブ出力 (``[output.photosphere]``) を読む (R2D2plus DEC-641、形式の正本は
R2D2plus ``docs/formats/photosphere.md``)。

光球近傍の鉛直範囲だけを高頻度に float32 で書く出力で、斜め視線のスペクトル合成や光球の磁場の解析に使う。
通常出力 (legacy / compressed) とも restart とも別の形式である。ここでは書かない。

配置::

    <directory>/format.toml               版・量の並び・単位・格子・ランクごとの範囲 (全スナップショット共通)
    <directory>/x.f64, y.f64, z.f64       物理セルの中心座標 (float64、cm。z は書いた鉛直範囲だけ)
    <directory>/<index:08d>/meta.toml     これがあれば書き終わっている (step・time)
    <directory>/<index:08d>/state.<rank:04d>.f32

各 state ファイルは、``format.toml`` の ``fields`` の順に量を並べ (field-major)、各量はランクの範囲を
(z, x, y) の C 順 (y が最速) で持つ little endian の float32。

**軸**: 配列は R2D2plus の内部の軸 ``(z, x, y)`` で、**z が鉛直**。この並びは pyR2D2 の legacy / compressed
出力 (Fortran の軸 ``(x, y, z)``、x が鉛直) と同じ順番である (R2D2plus の z・x・y が Fortran の x・y・z)。
量の名前の ``vertical``・``horizontal_x``・``horizontal_y`` は R2D2plus の z・x・y で、Fortran の名前では
``velocity_vertical`` = vx、``velocity_horizontal_x`` = vy、``velocity_horizontal_y`` = vz (磁場も同じ)。

版:

- 1: 6 量 (temperature, pressure, density, velocity_vertical, velocity_horizontal_x, velocity_horizontal_y) 固定。
- 2: ``fields`` に書いた量 (磁場 magnetic_vertical・magnetic_horizontal_x・magnetic_horizontal_y を含みうる)。
  単位は ``units`` (磁場は G)。
"""

from dataclasses import dataclass, field
from pathlib import Path
import tomllib

import numpy as np

__all__ = [
    "PHOTOSPHERE_FIELDS",
    "PHOTOSPHERE_UNITS",
    "Photosphere",
    "PhotosphereSnapshot",
    "read_photosphere",
]

#: 書きうる量の名前と、ファイル内の並びの基準 (正準順)。先頭 6 量が版 1。
PHOTOSPHERE_FIELDS = (
    "temperature",
    "pressure",
    "density",
    "velocity_vertical",
    "velocity_horizontal_x",
    "velocity_horizontal_y",
    "magnetic_vertical",
    "magnetic_horizontal_x",
    "magnetic_horizontal_y",
)
#: 量ごとの単位 (CGS、磁場は G)。版 1 の format.toml には ``units`` が無いのでこれを使う。
PHOTOSPHERE_UNITS = dict(zip(PHOTOSPHERE_FIELDS, ("K", "dyn/cm^2", "g/cm^3", "cm/s", "cm/s", "cm/s", "G", "G", "G")))

_SUPPORTED_VERSIONS = (1, 2)


@dataclass
class PhotosphereSnapshot:
    """:meth:`Photosphere.read` (:func:`read_photosphere`) の結果。"""

    path: Path
    index: int
    step: int
    time: float
    #: ``{name: ndarray (nz, nx, ny) float32}``。z が鉛直 (R2D2plus の軸)。
    fields: dict
    #: 読んだ範囲の座標 (cm、float64)
    x: np.ndarray
    y: np.ndarray
    z: np.ndarray
    units: dict = field(default_factory=dict)

    def __getitem__(self, name):
        return self.fields[name]


def _as_slice(value, n, axis):
    if value is None:
        return 0, n
    if isinstance(value, slice):
        start, stop, step = value.indices(n)
        if step != 1:
            raise ValueError(f"{axis}: only contiguous ranges (step 1) are supported")
        return start, max(start, stop)
    start, stop = value
    start, stop, _ = slice(start, stop).indices(n)
    return start, max(start, stop)


class Photosphere:
    """
    光球スラブ出力のディレクトリ。``format.toml`` と座標を一度だけ読む。

    Parameters
    ----------
    path : str or pathlib.Path
        ``format.toml`` のあるディレクトリ (既定の書き先は ``data/photosphere``)。ランのディレクトリや
        ``data/`` を渡すと、その下の ``photosphere`` を探す。

    Examples
    --------
    .. code-block:: python

        ph = pyR2D2.Photosphere("run/d003/data/photosphere")
        ph.fields, ph.units, ph.shape     # 書かれている量・単位・(nz, nx, ny)
        ph.snapshots()                    # 書き終わった番号
        s = ph.read(ph.snapshots()[-1], fields=["magnetic_vertical", "velocity_vertical"])
        s.time, s["magnetic_vertical"][-1]   # 最上段の鉛直磁場 (G)、shape (nx, ny)
    """

    def __init__(self, path):
        path = Path(path)
        for candidate in (path, path / "photosphere", path / "data" / "photosphere"):
            if (candidate / "format.toml").exists():
                path = candidate
                break
        else:
            raise FileNotFoundError(f"{path}: no format.toml of a photosphere output")
        self.path = path
        self.format = tomllib.loads((path / "format.toml").read_text())
        fmt = self.format["format"]
        self.version = int(fmt["version"])
        if self.version not in _SUPPORTED_VERSIONS:
            raise ValueError(f"{path}: unsupported photosphere format version {self.version} "
                             f"(this pyR2D2 reads {_SUPPORTED_VERSIONS})")
        if fmt.get("dtype") != "float32" or fmt.get("byte_order") != "little":
            raise ValueError(f"{path}: expected little endian float32, got {fmt.get('dtype')}/{fmt.get('byte_order')}")
        if list(fmt.get("axes", [])) != ["z", "x", "y"]:
            raise ValueError(f"{path}: expected axes [z, x, y], got {fmt.get('axes')}")
        #: 書かれている量の名前 (ファイル内の順)
        self.fields = tuple(fmt["fields"])
        unknown = [f for f in self.fields if f not in PHOTOSPHERE_FIELDS]
        if unknown:
            raise ValueError(f"{path}: unknown fields {unknown}")
        if self.version == 1 and self.fields != PHOTOSPHERE_FIELDS[:6]:
            raise ValueError(f"{path}: version 1 must have the 6 fixed fields, got {self.fields}")
        units = fmt.get("units")
        #: ``{name: 単位}``
        self.units = dict(zip(self.fields, units)) if units else {f: PHOTOSPHERE_UNITS[f] for f in self.fields}
        #: 全体の形 ``(nz, nx, ny)``
        self.shape = tuple(int(n) for n in self.format["grid"]["shape"])
        self.x = np.fromfile(path / self.format["grid"]["x"], dtype="<f8")
        self.y = np.fromfile(path / self.format["grid"]["y"], dtype="<f8")
        self.z = np.fromfile(path / self.format["grid"]["z"], dtype="<f8")
        if (self.z.size, self.x.size, self.y.size) != self.shape:
            raise ValueError(f"{path}: coordinates {(self.z.size, self.x.size, self.y.size)} do not match shape {self.shape}")
        self.ranks = [
            (int(r["id"]), tuple(int(n) for n in r["start"]), tuple(int(n) for n in r["shape"]))
            for r in self.format["rank"]
        ]
        self.rank_file_pattern = self.format["mpi"].get("rank_files", "state.{rank:04d}.f32")

    def snapshots(self):
        """書き終わった (``meta.toml`` がある) スナップショットの番号を昇順で返す。"""
        out = []
        for p in self.path.iterdir():
            if p.is_dir() and p.name.isdigit() and (p / "meta.toml").exists():
                out.append(int(p.name))
        return sorted(out)

    def read(self, index, fields=None, z=None, x=None, y=None):
        """
        1 つのスナップショットを読み、ランクのファイルを全体の配列に組み立てる。

        Parameters
        ----------
        index : int
            スナップショットの番号 (ディレクトリ名、``snapshots()`` の値)
        fields : list of str, optional
            読む量 (既定は書かれている全部)。書かれていない量を指定すると例外。
        z, x, y : slice or (start, stop), optional
            読む範囲 (スラブの中の添字、step 1 だけ)。``z=slice(-10, None)`` で上の 10 段など。
            必要なランクの必要な部分だけを読む (memmap)。

        Returns
        -------
        PhotosphereSnapshot
            ``fields[name]`` は shape ``(nz, nx, ny)`` の float32 (z が鉛直)。
        """
        names = self.fields if fields is None else tuple(fields)
        missing = [f for f in names if f not in self.fields]
        if missing:
            raise KeyError(f"{self.path}: fields {missing} are not in this output (have {self.fields})")
        snap = self.path / f"{int(index):08d}"
        if not (snap / "meta.toml").exists():
            raise FileNotFoundError(f"{snap}: no meta.toml (snapshot missing or incomplete)")
        meta = tomllib.loads((snap / "meta.toml").read_text())["snapshot"]
        if int(meta.get("format_version", 1)) != self.version:
            raise ValueError(f"{snap}: meta.toml format_version {meta.get('format_version')} "
                             f"differs from format.toml version {self.version}")

        nz, nx, ny = self.shape
        bounds = [_as_slice(z, nz, "z"), _as_slice(x, nx, "x"), _as_slice(y, ny, "y")]
        out_shape = tuple(b - a for a, b in bounds)
        out = {name: np.empty(out_shape, dtype=np.float32) for name in names}
        positions = [self.fields.index(name) for name in names]
        nfields = len(self.fields)

        for rank, start, shape in self.ranks:
            if 0 in shape:
                continue
            lo = [max(a, s) for (a, _), s in zip(bounds, start)]
            hi = [min(b, s + n) for (_, b), s, n in zip(bounds, start, shape)]
            if any(h <= l for l, h in zip(lo, hi)):
                continue
            file = snap / self.rank_file_pattern.format(rank=rank)
            expected = nfields * shape[0] * shape[1] * shape[2] * 4
            if file.stat().st_size != expected:
                raise ValueError(f"{file}: {file.stat().st_size} bytes, expected {expected} "
                                 f"({nfields} fields x {shape})")
            data = np.memmap(file, dtype="<f4", mode="r", shape=(nfields, *shape))
            src = tuple(slice(l - s, h - s) for l, h, s in zip(lo, hi, start))
            dst = tuple(slice(l - a, h - a) for l, h, (a, _) in zip(lo, hi, bounds))
            for name, n in zip(names, positions):
                out[name][dst] = data[(n, *src)]
            del data

        return PhotosphereSnapshot(
            path=snap,
            index=int(meta.get("index", index)),
            step=int(meta["step"]),
            time=float(meta["time"]),
            fields=out,
            x=self.x[bounds[1][0]:bounds[1][1]].copy(),
            y=self.y[bounds[2][0]:bounds[2][1]].copy(),
            z=self.z[bounds[0][0]:bounds[0][1]].copy(),
            units={name: self.units[name] for name in names},
        )


def read_photosphere(path, index="latest", fields=None, z=None, x=None, y=None):
    """
    光球スラブ出力の 1 スナップショットを読む (:class:`Photosphere` の ``read`` の近道)。

    Parameters
    ----------
    path : str or pathlib.Path
        出力のディレクトリ (``format.toml`` のあるところ)、``data/``、またはランのディレクトリ
    index : int or "latest"
        スナップショットの番号。``"latest"`` (既定) は書き終わった最後の番号。
    fields, z, x, y
        :meth:`Photosphere.read` と同じ。

    Returns
    -------
    PhotosphereSnapshot
    """
    ph = Photosphere(path)
    if index == "latest":
        have = ph.snapshots()
        if not have:
            raise FileNotFoundError(f"{ph.path}: no complete snapshot")
        index = have[-1]
    return ph.read(index, fields=fields, z=z, x=x, y=y)
