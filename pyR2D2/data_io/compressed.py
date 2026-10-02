"""R2D2plus の新出力形式 ``format = "compressed"`` (R2D2plus DEC-576) の読み込み。

形式の見分け方
--------------
``<datadir>/param/format.toml`` があり ``format = "compressed"`` なら新形式、
無ければ従来形式 (Fortran R2D2 と R2D2plus の ``format = "legacy"``) である。
:py:class:`pyR2D2.Data` はこれを自動で見分けるので、利用者は従来どおり
``d.qf.read(n)`` などを呼べばよい。

新形式の要点
------------
* 3D の ``remap/qq`` は 1 ランク 1 ファイルの ``r2d2plus-z`` コンテナ
  (変数ごとに byte shuffle + zstd)。中身の配列は従来形式とビット単位で同一。
* ``pr``, ``te``, ``op`` (と既定では ``ph``) は書かれない。pyR2D2 は要求された
  ときに ``ro``, ``se`` と EOS 表 (``param/eos_table.tbl``) から、R2D2plus の
  従来形式の書き出しと同じ式・同じ精度 (単精度の丸めの位置まで) で計算する。
  式は出力の種類ごとに違う (:class:`CompressedEOS` を参照)。
* tau / slice / 2D remap は従来と同じ単精度の生バイナリで、変数の並びだけが違う。
* prev / aftr は単精度の ``r2d2plus-z`` (従来は倍精度の生バイナリ)。

このモジュールの公開関数
------------------------
* :func:`read_r2d2plus_z` : ``r2d2plus-z`` ファイル 1 個を読む
* :func:`read_format_toml` : ``param/format.toml`` を読む
* :func:`write_initial_state` : R2D2plus の初期状態ファイル (``[initial_condition] type = "file"``) を書く
* :func:`write_boundary_series` : R2D2plus の境界条件 ``"file"`` 用の境界値の時系列を書く
* :func:`read_restart_meta`, :func:`read_restart_state` : R2D2plus の
  リスタート (``restart/<slot>``) を読む (**実験的**)
"""

import os
from pathlib import Path

import numpy as np

try:  # Python 3.11+
    import tomllib
except ImportError:  # pragma: no cover
    tomllib = None

R2D2PLUS_Z_MAGIC = b"R2D2PLUS-Z 1\n"

# format.toml が無い項目の既定値 (R2D2plus DEC-576 の version 1)
_FORMAT_DEFAULTS = {
    "variables": ["ro", "vx", "vy", "vz", "bx", "by", "bz", "se"],
    "eos_table": "eos_table.tbl",
    "remap_qq": "remap/qq/{group:05d}/{rank:08d}/qq.z.{nd:08d}.{rank:08d}",
    "tau": "tau/qq.c.{nd_tau:08d}",
    "tau_quantities": [
        "in", "ro", "se", "vx", "vy", "vz", "bx", "by", "bz", "height", "fr",
    ],
    "tau_height_origin": "rstar",
    "slice": "slice/{name}.c.{nd_tau:08d}.{n:08d}",
    "remap_2d": "remap/qq/qq.c.{nd:08d}",
    "remap_2d_extra": ["tu", "fr"],
    "prev": "prev/{group:05d}/{rank:08d}/qq.z.{nd:08d}.{index:08d}.{rank:08d}",
    "aftr": "aftr/{group:05d}/{rank:08d}/qq.z.{nd:08d}.{index:08d}.{rank:08d}",
}

# ランクの小ディレクトリは 1000 ランクごと (R2D2 と同じ)
RANK_GROUP_SIZE = 1000


def _load_toml(text):
    if tomllib is None:  # pragma: no cover
        raise ImportError("reading R2D2plus compressed output needs Python >= 3.11 (tomllib)")
    return tomllib.loads(text)


def read_format_toml(datadir):
    """
    ``<datadir>/param/format.toml`` を読む。

    Parameters
    ----------
    datadir : str or pathlib.Path
        ラン出力のディレクトリ (``.../data``)

    Returns
    -------
    info : dict or None
        ファイルが無ければ ``None`` (従来形式)。あれば中身の dict に、
        書かれていない項目の既定値を補ったもの。
    """
    path = Path(datadir) / "param" / "format.toml"
    if not path.exists():
        return None
    with open(path, "rb") as f:
        info = _load_toml(f.read().decode("utf-8"))
    for key, value in _FORMAT_DEFAULTS.items():
        info.setdefault(key, value)
    return info


# ---------------------------------------------------------------------------
# r2d2plus-z コンテナ
# ---------------------------------------------------------------------------


def _zstd_decompress(data, raw_bytes):
    """zstd を解く。使えるものを順に試す (Python 3.14 の標準、numcodecs、zstandard)。"""
    try:
        from compression import zstd as _std_zstd  # Python 3.14+

        out = _std_zstd.decompress(data)
    except ImportError:
        try:
            import numcodecs

            out = numcodecs.Zstd().decode(data)
        except ImportError:
            try:
                import zstandard
            except ImportError:
                raise ImportError(
                    "reading R2D2plus compressed output (r2d2plus-z) needs a zstd "
                    "decoder: install numcodecs (`pip install numcodecs`, also "
                    "installed with pyR2D2[zarr]) or zstandard"
                ) from None
            out = zstandard.ZstdDecompressor().decompress(data, max_output_size=raw_bytes)
    out = memoryview(out).cast("B")
    if out.nbytes != raw_bytes:
        raise ValueError(f"r2d2plus-z: decoded {out.nbytes} bytes, header says {raw_bytes}")
    return out


def _unshuffle(shuffled, dtype):
    """byte shuffle を戻す (要素のバイト位置ごとの面 → 要素の並び)。"""
    dtype = np.dtype(dtype)
    planes = np.frombuffer(shuffled, dtype=np.uint8).reshape(dtype.itemsize, -1)
    return planes.T.copy().view(dtype).reshape(-1)


def read_r2d2plus_z_header(path):
    """
    ``r2d2plus-z`` ファイルのヘッダだけを読む。

    Returns
    -------
    header : dict
        TOML ヘッダの中身 (``kind``, ``attributes``, ``variables`` など)。
        ``header["data_offset"]`` にデータ部の先頭のバイト位置を足してある。
    """
    with open(path, "rb") as f:
        magic = f.readline()
        if magic != R2D2PLUS_Z_MAGIC:
            raise ValueError(f"{path}: not a r2d2plus-z file")
        length = int(f.readline().strip(), 16)
        text = f.read(length)
        if len(text) != length:
            raise ValueError(f"{path}: truncated r2d2plus-z header")
        data_offset = f.tell()
    header = _load_toml(text.decode("utf-8"))
    if header.get("format") != "r2d2plus-z" or header.get("version") != 1:
        raise ValueError(f"{path}: unsupported r2d2plus-z format/version")
    header.setdefault("attributes", {})
    header.setdefault("variables", [])
    header["data_offset"] = data_offset
    return header


def _xxh3_hex(raw):
    try:
        import xxhash
    except ImportError:
        return None
    return f"{xxhash.xxh3_64_intdigest(raw):016x}"


def read_r2d2plus_z(path, names=None, verify=None):
    """
    R2D2plus の ``r2d2plus-z`` ファイル (新出力形式の remap/qq, prev, aftr、
    圧縮リスタート) を読む。

    Parameters
    ----------
    path : str or pathlib.Path
        ファイル
    names : iterable of str, optional
        読む変数名。省略時は全部。
    verify : bool or None, optional
        各変数の XXH3-64 を照合するか。``None`` (既定) なら ``xxhash``
        モジュールがあるときだけ照合する。``True`` で ``xxhash`` が無ければ例外。

    Returns
    -------
    header : dict
        TOML ヘッダ (:func:`read_r2d2plus_z_header` と同じ)
    arrays : dict of numpy.ndarray
        ``{name: array}``。ヘッダの ``shape`` と ``order`` どおりの形で、
        dtype はヘッダの ``dtype`` (little endian)。

    Examples
    --------
    .. code-block:: python

        header, qq = pyR2D2.data_io.compressed.read_r2d2plus_z(
            "data/remap/qq/00000/00000000/qq.z.00000001.00000000")
        qq["ro"].shape   # (iixl, jjxl, kx), 鉛直が最速 (Fortran 順)
    """
    header = read_r2d2plus_z_header(path)
    variables = header["variables"]
    if names is not None:
        names = set(names)
        unknown = names - {v["name"] for v in variables}
        if unknown:
            raise KeyError(f"{path}: no variable(s) {sorted(unknown)}")
        variables = [v for v in variables if v["name"] in names]

    if verify and _xxh3_hex(b"") is None:
        raise ImportError("verify=True needs the xxhash module (`pip install xxhash`)")

    arrays = {}
    with open(path, "rb") as f:
        for v in variables:
            if v.get("codec") != "shuffle+zstd":
                raise ValueError(f"{path}: unsupported codec {v.get('codec')!r}")
            dtype = {"float32": "<f4", "float64": "<f8"}.get(v["dtype"])
            if dtype is None:
                raise ValueError(f"{path}: unsupported dtype {v['dtype']!r}")
            f.seek(header["data_offset"] + int(v["offset"]))
            stored = f.read(int(v["stored_bytes"]))
            if len(stored) != int(v["stored_bytes"]):
                raise ValueError(f"{path}: truncated data for {v['name']}")
            shuffled = _zstd_decompress(stored, int(v["raw_bytes"]))
            flat = _unshuffle(shuffled, dtype)
            if verify is not False:
                digest = _xxh3_hex(flat.tobytes())
                if digest is not None and digest != v["xxh3"]:
                    raise ValueError(
                        f"{path}: XXH3 mismatch for {v['name']} ({digest} != {v['xxh3']})"
                    )
            arrays[v["name"]] = flat.reshape(tuple(v["shape"]), order=v.get("order", "C"))
    return header, arrays


def _zstd_compress(data, level):
    """zstd で圧縮する (読み口と同じく、使えるものを順に試す)。"""
    try:
        from compression import zstd as _std_zstd  # Python 3.14+

        return _std_zstd.compress(data, level=level)
    except ImportError:
        pass
    try:
        import numcodecs

        return bytes(numcodecs.Zstd(level=level).encode(data))
    except ImportError:
        pass
    try:
        import zstandard
    except ImportError:
        raise ImportError(
            "writing a r2d2plus-z file needs a zstd encoder: install numcodecs "
            "(`pip install numcodecs`) or zstandard"
        ) from None
    return zstandard.ZstdCompressor(level=level).compress(data)


def _require_xxh3():
    try:
        import xxhash
    except ImportError:
        raise ImportError(
            "writing a r2d2plus-z file needs the xxhash module (`pip install xxhash`, "
            "or pyR2D2[compressed]): R2D2plus verifies the XXH3-64 checksum of every "
            "variable, and there is no practical pure-Python XXH3"
        ) from None
    return xxhash


def _toml_scalar(key, value):
    """[attributes] の 1 行。R2D2plus が受け付けるのはスカラー (文字列・整数・実数・真偽)。"""
    import json
    import math

    if not isinstance(key, str) or not key or not all(c.isalnum() or c == "_" for c in key) \
            or not key.isascii():
        raise ValueError(f"attribute name must be [A-Za-z0-9_]+, got {key!r}")
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, bool):
        text = "true" if value else "false"
    elif isinstance(value, int):
        text = str(value)
    elif isinstance(value, float):
        if not math.isfinite(value):
            raise ValueError(f"attribute {key!r} must be finite, got {value}")
        text = repr(value)
        if "." not in text and "e" not in text and "E" not in text:
            text += ".0"
    elif isinstance(value, (str, os.PathLike)):
        text = json.dumps(os.fspath(value))  # JSON の文字列は TOML の basic string でもある
    else:
        raise TypeError(f"attribute {key!r} must be a str/int/float/bool scalar, got {type(value)}")
    return f"{key} = {text}\n"


def write_r2d2plus_z(path, kind, variables, attributes=None, level=3):
    """
    ``r2d2plus-z`` ファイルを書く (R2D2plus ``CompressedFileWriter`` と同じ形)。

    Parameters
    ----------
    path : str or pathlib.Path
        書き出し先。一時ファイルに書いてから置き換える。
    kind : str
        ヘッダの ``kind``
    variables : list of (name, numpy.ndarray, order)
        変数。dtype は float32 か float64 (little endian で書く)、``order`` は
        ``"F"`` か ``"C"`` (``shape`` はその並びで解釈される)。
    attributes : dict, optional
        ``[attributes]`` に書くスカラー
    level : int
        zstd の圧縮レベル

    Notes
    -----
    各変数の XXH3-64 (R2D2plus が読むときに照合する) の計算に ``xxhash`` が要る。
    """
    xxhash = _require_xxh3()
    blocks, entries, offset = [], [], 0
    seen = set()
    for name, array, order in variables:
        if not isinstance(name, str) or not name or len(name) > 64 or not name.isascii() \
                or not all(c.isalnum() or c == "_" for c in name) or name in seen:
            raise ValueError(f"bad or duplicate variable name {name!r}")
        seen.add(name)
        if order not in ("F", "C"):
            raise ValueError(f"order must be 'F' or 'C', got {order!r}")
        array = np.asarray(array)
        dtype_name = {"f4": "float32", "f8": "float64"}.get(array.dtype.str[1:])
        if array.dtype.kind != "f" or dtype_name is None:
            raise TypeError(f"{name}: dtype must be float32 or float64, got {array.dtype}")
        raw = np.asarray(array, dtype=array.dtype.newbyteorder("<")).tobytes(order=order)
        es = array.dtype.itemsize
        shuffled = np.frombuffer(raw, np.uint8).reshape(-1, es).T.tobytes()
        stored = _zstd_compress(shuffled, level)
        blocks.append(stored)
        entries.append(
            "\n[[variables]]\n"
            f'name = "{name}"\ndtype = "{dtype_name}"\norder = "{order}"\n'
            f"shape = [{', '.join(str(int(s)) for s in array.shape)}]\n"
            'codec = "shuffle+zstd"\n'
            f"offset = {offset}\nstored_bytes = {len(stored)}\nraw_bytes = {len(raw)}\n"
            f'xxh3 = "{xxhash.xxh3_64_intdigest(raw):016x}"\n'
        )
        offset += len(stored)
    header = f'format = "r2d2plus-z"\nversion = 1\nkind = "{kind}"\n\n[attributes]\n'
    header += "".join(_toml_scalar(k, v) for k, v in (attributes or {}).items())
    header += "".join(entries)
    header = header.encode("utf-8")

    path = Path(path)
    temporary = path.with_name(path.name + ".tmp")
    with open(temporary, "wb") as f:
        f.write(R2D2PLUS_Z_MAGIC)
        f.write(f"{len(header):016x}\n".encode())
        f.write(header)
        for b in blocks:
            f.write(b)
    os.replace(temporary, path)
    return path


INITIAL_STATE_VARIABLES = ["ro", "vx", "vy", "vz", "bx", "by", "bz", "se", "ph"]


def write_initial_state(path, fields, attributes=None, level=3):
    """
    R2D2plus の初期状態ファイル (``r2d2plus-z``、``kind = "initial_state"``) を書く。

    R2D2plus では次のように指定して読む::

        [initial_condition]
        type = "file"
        [initial_condition.params]
        path = "/path/to/initial_state.z"

    **ro と se は背景からの摂動である。** 従来形式の ``d.qf.ro`` / ``d.qf.se``
    (``remap/qq``) が返すものと同じ量で、全量 (``ro0 + ro``) を渡してはならない。
    速度・磁場・``ph`` はそのままの値。

    Parameters
    ----------
    path : str or pathlib.Path
        書き出し先
    fields : dict of numpy.ndarray
        ``{name: array}``。名前は ``ro, vx, vy, vz, bx, by, bz, se, ph`` から。
        **名前と軸は pyR2D2 の 3D 配列 (Fortran の並び) と同じ**: ``vx``, ``bx`` が
        鉛直成分、配列の形は ``(ix, jx, kx)`` = (鉛直, 第1水平, 第2水平)、
        領域全体の物理セルのみ (margin 無し)。``d.qf.read(n)`` で読んだ配列を
        そのまま渡せばよい。すべて同じ形であること。
    attributes : dict, optional
        ``[attributes]`` に書くスカラー (文字列・整数・実数・真偽)。例:
        ``{"source": "d001", "nd": 120, "time": 3600.0}``
    level : int, optional
        zstd の圧縮レベル (既定 3)

    Returns
    -------
    path : pathlib.Path

    Notes
    -----
    各変数は float64、``order = "F"`` (鉛直が最速)、``shape = list(array.shape)``
    で書く。XXH3-64 のチェックサムを R2D2plus が照合するので ``xxhash`` が必要
    (``pip install xxhash``)。zstd には numcodecs (または zstandard) を使う。

    Examples
    --------
    前の出力を初期条件にする (従来形式・新形式どちらのランでも同じ)::

        import pyR2D2
        d = pyR2D2.Data("../run/d001/data")
        n = 120
        d.qf.read(n, keys=["ro", "vx", "vy", "vz", "bx", "by", "bz", "se"])
        fields = {k: d.qf.__dict__[k] for k in
                  ["ro", "vx", "vy", "vz", "bx", "by", "bz", "se"]}
        pyR2D2.write_initial_state(
            "initial_state.z", fields,
            attributes={"source": str(d.datadir), "nd": n,
                        "time": float(d.time_read(n, verbose=False))})
    """
    if not fields:
        raise ValueError("fields is empty")
    unknown = [k for k in fields if k not in INITIAL_STATE_VARIABLES]
    if unknown:
        raise ValueError(f"unknown variable(s) {unknown}; allowed: {INITIAL_STATE_VARIABLES}")
    shape = None
    variables = []
    for name in INITIAL_STATE_VARIABLES:  # 並びは固定 (Fortran の成分順)
        if name not in fields:
            continue
        array = np.asarray(fields[name], dtype=np.float64)
        if array.ndim != 3:
            raise ValueError(f"{name}: expected a 3D array (ix, jx, kx), got shape {array.shape}")
        if shape is None:
            shape = array.shape
        elif array.shape != shape:
            raise ValueError(f"{name}: shape {array.shape} differs from {shape}")
        if not np.all(np.isfinite(array)):
            raise ValueError(f"{name}: contains NaN or Inf")
        variables.append((name, array, "F"))
    if "ro" in fields and np.all(variables[0][1] > 0) and variables[0][0] == "ro":
        import warnings

        warnings.warn(
            "ro is positive everywhere: initial-state ro must be the perturbation from "
            "the background (as d.qf.ro), not the total density",
            stacklevel=2,
        )
    return write_r2d2plus_z(path, "initial_state", variables, attributes, level)


def write_boundary_series(path, times, values, attributes=None, level=3):
    """
    R2D2plus の境界値の時系列ファイル (``r2d2plus-z``、``kind = "boundary_series"``) を書く。

    R2D2plus の汎用境界条件の規則 ``"file"`` が読む。例::

        [boundary.params]
        bottom.vz = "file"
        bottom.vz_file = "/path/to/bottom_vz.z"

    値はその壁の ghost 層すべてに使われ、時間方向にはステップ開始時刻で線形内挿
    される (範囲外は端の値で止める)。

    **[boundary.params] の変数名は C++ の成分名である** (``ro, vx, vy, vz, bx, by,
    bz, se, ps``、**C++ では z が鉛直**)。pyR2D2 の 3D 配列 (Fortran の名前) とは
    速度と磁場の名前がずれる:

    ============================  ==========================
    pyR2D2 / Fortran の名前       R2D2plus [boundary.params]
    ============================  ==========================
    ``vx`` (鉛直)                 ``vz``
    ``vy`` (第1水平)              ``vx``
    ``vz`` (第2水平)              ``vy``
    ``bx``, ``by``, ``bz``        ``bz``, ``bx``, ``by``
    ``ph``                        ``ps``
    ``ro``, ``se``                ``ro``, ``se`` (同じ)
    ============================  ==========================

    つまり鉛直速度を与えるなら ``d.qf.vx`` 由来の値を ``bottom.vz`` に指定する。
    **ro と se の値は背景からの摂動** (``d.qf.ro``, ``d.qf.se`` と同じ量)。

    Parameters
    ----------
    path : str or pathlib.Path
        書き出し先
    times : array_like, shape (nt,)
        時刻 [s]。狭義単調増加であること。
    values : numpy.ndarray, shape (n1, n2, nt)
        境界値。``n1`` = 第1水平 (pyR2D2 の ``jx``、Fortran y = C++ x)、
        ``n2`` = 第2水平 (``kx``、Fortran z = C++ y)、領域全体の物理セルのみ。
        **時刻は最後の軸**。例えば ``d.qf.vx[0]`` (下端の面、形 ``(jx, kx)``) を
        時刻ごとに集めて ``np.stack(planes, axis=-1)`` とすればよい。
    attributes : dict, optional
        ``[attributes]`` に書くスカラー
    level : int, optional
        zstd の圧縮レベル (既定 3)

    Returns
    -------
    path : pathlib.Path

    Notes
    -----
    変数は ``time`` (float64, ``[nt]``) と ``value`` (float64, order ``"F"``、
    ``[n1, n2, nt]``、第1水平が最速) の 2 つ。XXH3-64 の計算に ``xxhash`` が要る。
    """
    times = np.asarray(times, dtype=np.float64)
    values = np.asarray(values, dtype=np.float64)
    if times.ndim != 1 or times.size == 0:
        raise ValueError(f"times must be a non-empty 1D array, got shape {times.shape}")
    if not np.all(np.isfinite(times)) or np.any(np.diff(times) <= 0):
        raise ValueError("times must be finite and strictly increasing")
    if values.ndim != 3 or values.shape[2] != times.size:
        raise ValueError(
            f"values must have shape (n1, n2, nt) with nt = {times.size} (time last), "
            f"got {values.shape}"
        )
    if not np.all(np.isfinite(values)):
        raise ValueError("values contains NaN or Inf")
    return write_r2d2plus_z(
        path,
        "boundary_series",
        [("time", times, "F"), ("value", values, "F")],
        attributes,
        level,
    )


# ---------------------------------------------------------------------------
# 派生量 (pr, te, op) を R2D2plus の従来形式の書き出しと同じ式で作る
# ---------------------------------------------------------------------------


class CompressedEOS:
    """
    新形式が書かない ``pr``, ``te``, ``op`` を ``ro``, ``se`` から計算する。

    R2D2plus の従来形式の書き出しは出力の種類ごとに式が違うので、それぞれを
    **同じ式・同じ精度・同じ丸めの位置**で写してある。表の内挿は R2D2plus
    ``EosTableView::interpolate`` (系統 R) と同じ: 添字は除算で求めて 0 方向へ
    切り捨て、両端でクランプ、``log(値)`` を (log ρ, s) で双線形内挿し、総和の
    あとに刻みの逆数を 2 回掛けて ``exp`` する。対数は ``log(x + 1e-200)``、
    不透明度だけ ``log(max(x, 1e-20))``。

    ``pyR2D2.util.eos_table`` / ``eos_switch`` とは総和のスケールの仕方
    (割るか逆数を掛けるか) と不透明度の対数の取り方が違い、最終桁が変わる
    ことがあるので、ここでは使わない。

    Parameters
    ----------
    table : dict
        EOS 表の生の値 (``ro`` = log ρ 軸、``se`` 軸、``pr``, ``te``, ``op``
        の 2 次元配列 (n_ro, n_se))
    param : pyR2D2.Parameters
        背景量 (``ro0``, ``se0``, ``pr0``, ``te0``, ``dprdro``, ``dprdse``,
        ``dtedro``, ``dtedse``) を持つもの
    """

    quantities = ("pr", "te", "op")

    def __init__(self, table, param):
        self.log_rho_axis = np.asarray(table["ro"], dtype=np.float64)
        self.se_axis = np.asarray(table["se"], dtype=np.float64)
        self.d_log_rho = self.log_rho_axis[1] - self.log_rho_axis[0]
        self.d_se = self.se_axis[1] - self.se_axis[0]
        self.d_log_rho_inv = 1.0 / self.d_log_rho
        self.d_se_inv = 1.0 / self.d_se
        self.log_value = {}
        for q in self.quantities:
            v = np.asarray(table[q], dtype=np.float64)
            if q == "op":
                self.log_value[q] = np.log(np.maximum(v, 1.0e-20))
            else:
                self.log_value[q] = np.log(v + 1.0e-200)

        # 背景 (margin を除いた鉛直 1D)。倍精度と、単精度へ落としたもの
        # (3D remap の書き出しは単精度の背景を使う)。
        names = ["ro0", "se0", "pr0", "te0", "dprdro", "dprdse", "dtedro", "dtedse"]
        self.bg64 = {k: np.asarray(getattr(param, k), dtype=np.float64) for k in names}
        self.bg32 = {k: v.astype(np.float32) for k, v in self.bg64.items()}

    def interpolate(self, q, rho, se):
        """
        表を引く (倍精度の総量 ``rho``, ``se`` から、倍精度の値を返す)。
        R2D2plus ``EosAxes::locate`` + ``EosTableView::interpolate`` と同じ演算順。
        """
        rho = np.asarray(rho, dtype=np.float64)
        se = np.asarray(se, dtype=np.float64)
        with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
            log_rho = np.log(rho)
            fi = (log_rho - self.log_rho_axis[0]) / self.d_log_rho
            fj = (se - self.se_axis[0]) / self.d_se
            # C++ の static_cast<int> (0 方向への切り捨て) のあとクランプ
            i = np.clip(np.nan_to_num(np.trunc(fi)), 0, self.log_rho_axis.size - 2).astype(np.intp)
            j = np.clip(np.nan_to_num(np.trunc(fj)), 0, self.se_axis.size - 2).astype(np.intp)
            dro0 = log_rho - self.log_rho_axis[i]
            dro1 = self.d_log_rho - dro0
            dse0 = se - self.se_axis[j]
            dse1 = self.d_se - dse0
            v = self.log_value[q]
            # 足し算の順と括弧の位置は C++ と同じ (最終桁が変わるので変えない)
            s = v[i, j] * dro1 * dse1
            s = s + v[i + 1, j] * dro0 * dse1
            s = s + v[i, j + 1] * dro1 * dse0
            s = s + v[i + 1, j + 1] * dro0 * dse0
            return np.exp(s * self.d_log_rho_inv * self.d_se_inv)

    # --- 3D remap/qq (R2D2plus legacy_state_writer.cpp) -----------------------
    def remap(self, q, ro, se, iz):
        """
        3D ``remap/qq`` の ``pr``, ``te`` (摂動), ``op`` (全量)。

        すべて単精度: ``ct = max(|ro|/ro0, |se|/se0)``、``feos`` は
        ``ct >= 3e-3`` で 1 (硬い切り替え)、線形側 ``dXdro*ro + dXdse*se``、
        表側 ``float32(interp(double(ro+ro0), double(se+se0)) - X0)``、
        ``X = table*feos + linear*(1-feos)``。``op`` は切り替え無しの表の値。

        Parameters
        ----------
        q : str
            ``"pr"``, ``"te"`` または ``"op"``
        ro, se : numpy.ndarray (float32)
            出力されている摂動
        iz : numpy.ndarray of int
            各点の鉛直の添字 (margin を除いた 0 始まり、``ro`` へ broadcast できる形)
        """
        f32 = np.float32
        bg = self.bg32
        ro = np.asarray(ro, dtype=f32)
        se = np.asarray(se, dtype=f32)
        ro0 = bg["ro0"][iz]
        se0 = bg["se0"][iz]
        rho_total = (ro + ro0).astype(np.float64)
        se_total = (se + se0).astype(np.float64)
        if q == "op":
            return self.interpolate("op", rho_total, se_total).astype(f32)
        with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
            ct = np.maximum(np.abs(ro) / ro0, np.abs(se) / se0)
            feos = (ct.astype(np.float64) - 3.0e-3 >= 0.0).astype(f32)
            linear = bg[f"d{q}dro"][iz] * ro + bg[f"d{q}dse"][iz] * se
            table = (
                self.interpolate(q, rho_total, se_total) - bg[f"{q}0"][iz].astype(np.float64)
            ).astype(f32)
            return table * feos + linear * (f32(1.0) - feos)

    # --- tau (R2D2plus tau_surface_writer.cpp) ------------------------------
    def total(self, q, ro_total, se_total):
        """
        tau 面の ``pr``, ``te`` (全量)。tau の ``ro``, ``se`` は全量の単精度で、
        R2D2plus はその単精度値を倍精度へ広げて表を引き、単精度へ丸めている
        (切り替え無し)。面が見つからなかった点 (``ro == se == 0``) は 0 にする。
        """
        ro_total = np.asarray(ro_total)
        se_total = np.asarray(se_total)
        value = self.interpolate(q, ro_total.astype(np.float64), se_total.astype(np.float64))
        value = value.astype(np.float32)
        value[(ro_total == 0) & (se_total == 0)] = 0.0
        return value

    # --- slice (R2D2plus slice_writer.cpp) ----------------------------------
    def slice(self, q, ro, se, iz):
        """
        slice の ``pr``, ``te`` (摂動)。R2D2plus は**倍精度**の状態から切り替え付き
        で作って単精度へ丸めている。ここでは単精度で出力された ``ro``, ``se``
        から同じ倍精度の式で作るので、入力の丸めの分だけ違いうる。
        """
        bg = self.bg64
        ro = np.asarray(ro, dtype=np.float64)
        se = np.asarray(se, dtype=np.float64)
        ro0 = bg["ro0"][iz]
        se0 = bg["se0"][iz]
        with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
            linear = bg[f"d{q}dro"][iz] * ro + bg[f"d{q}dse"][iz] * se
            ct = np.maximum(np.abs(ro) / ro0, np.abs(se) / se0)
            feos = np.where(ct - 3.0e-3 >= 0.0, 1.0, 0.0)
            table = self.interpolate(q, ro + ro0, se + se0) - bg[f"{q}0"][iz]
            return (table * feos + linear * (1.0 - feos)).astype(np.float32)

    # --- 2D remap (R2D2plus remap2d_writer.cpp) -----------------------------
    def remap2d(self, q, ro, se, iz):
        """
        2D ``remap/qq`` の ``pr``, ``te`` (摂動), ``op`` (全量)。R2D2plus の従来形式は
        **切り替え無しの表の値**を倍精度の状態から作っている。ここでは単精度で
        出力された ``ro``, ``se`` から同じ式で作る。
        """
        bg = self.bg64
        ro = np.asarray(ro, dtype=np.float64)
        se = np.asarray(se, dtype=np.float64)
        value = self.interpolate(q, ro + bg["ro0"][iz], se + bg["se0"][iz])
        if q != "op":
            value = value - bg[f"{q}0"][iz]
        return value.astype(np.float32)


def load_eos_table_raw(path):
    """EOS 表 (r2d2plus-table の ``.tbl`` か npz) を生の値の dict で読む。"""
    from .parameters import _read_r2d2plus_table

    with open(path, "rb") as f:
        magic = f.read(15)
    if magic == b"r2d2plus-table ":
        return _read_r2d2plus_table(path)
    with np.load(path) as npz:
        return {k: npz[k] for k in npz.files}


def compressed_eos(param):
    """``param`` (pyR2D2.Parameters) に対応する :class:`CompressedEOS` を作って覚える。"""
    eos = param.__dict__.get("_compressed_eos")
    if eos is not None:
        return eos
    path = param.__dict__.get("eos_table_path")
    if path is None or not os.path.exists(path):
        raise FileNotFoundError(
            "this run uses the R2D2plus compressed output format, which does not store "
            "pr/te/op; computing them needs the EoS table (expected at "
            f"{Path(param.datadir) / 'param' / param.format_info['eos_table']}). "
            "Read ro/se etc. only, or give the table with d.p.read_eos_table(path)."
        )
    eos = CompressedEOS(load_eos_table_raw(path), param)
    param.__dict__["_compressed_eos"] = eos
    return eos


def format_path(datadir, info, kind, **fields):
    """``format.toml`` のパスの雛形 ``info[kind]`` を埋めて ``datadir`` からのパスにする。"""
    if "rank" in fields:
        fields.setdefault("group", fields["rank"] // RANK_GROUP_SIZE)
    return Path(datadir) / info[kind].format(**fields)


# ---------------------------------------------------------------------------
# リスタート (実験的)
# ---------------------------------------------------------------------------

RESTART_COMPONENTS = ["ro", "vx", "vy", "vz", "bx", "by", "bz", "se", "ps"]


def read_restart_meta(slot_dir):
    """
    R2D2plus のリスタートのスロット (``data/restart/e`` など) の ``meta.toml`` を読む。

    **実験的**: リスタートは pyR2D2 の読み込み対象として設計された形式ではない
    (R2D2plus DEC-076)。形式が変わりうる。
    """
    with open(Path(slot_dir) / "meta.toml", "rb") as f:
        return _load_toml(f.read().decode("utf-8"))


def read_restart_state(slot_dir, rank, state_subdir="", verify=None):
    """
    R2D2plus のリスタートから 1 ランク分の状態を読む (**実験的**)。

    ``meta.toml`` の ``generation`` の下の ``state.{rank:08d}.bin`` を読む。
    中身は生の float64 か、``r2d2plus-z`` の変数 ``state`` (圧縮リスタート)。

    Parameters
    ----------
    slot_dir : str or pathlib.Path
        スロットのディレクトリ (``data/restart/e``, ``data/restart/o``,
        ``data/restart/00000000`` など)
    rank : int
        MPI ランク
    state_subdir : str, optional
        YinYang のパネルの小ディレクトリ (R2D2plus DEC-430)
    verify : bool or None, optional
        ``meta.toml`` の XXH3-64 (ファイル全体) を照合するか。``None`` なら
        ``xxhash`` があるときだけ。

    Returns
    -------
    state : dict of numpy.ndarray
        ``ro, vx, vy, vz, bx, by, bz, se, ps`` の 9 成分、各 ``(nx, ny, nz)``
        (そのランクの物理セルのみ、float64)。**軸と成分は C++ の並び**で、
        ``z`` (最後の軸) が鉛直、``vz``/``bz`` が鉛直成分である
        (Fortran の ``vx`` = C++ の ``vz``、``vy`` = ``vx``、``vz`` = ``vy``)。
    meta : dict
        ``meta.toml`` の中身
    """
    slot_dir = Path(slot_dir)
    meta = read_restart_meta(slot_dir)
    restart = meta.get("restart", {})
    grid = meta["grid"]
    mpi = meta["mpi"]
    base = slot_dir / restart["generation"] if restart.get("generation") else slot_dir
    if state_subdir:
        base = base / state_subdir
    path = base / f"state.{rank:08d}.bin"
    with open(path, "rb") as f:
        raw = f.read()

    checksums = meta.get("checksum", {}).get("values")
    # checksum は world rank 順 (YinYang は Yin 0..P-1、Yang P..2P-1。R2D2plus restart.hpp)。
    # 2026-10-02 まで Yang 側も rank で引いており、Yin 側の値と比べていた。
    world = rank + (mpi["size_x"] * mpi["size_y"] * mpi["size_z"] if state_subdir == "yang" else 0)
    if verify is not False and checksums is not None and world < len(checksums):
        digest = _xxh3_hex(raw)
        if digest is None and verify:
            raise ImportError("verify=True needs the xxhash module")
        if digest is not None and digest != checksums[world]:
            raise ValueError(f"{path}: XXH3 checksum mismatch ({digest} != {checksums[world]})")

    if raw.startswith(R2D2PLUS_Z_MAGIC):
        _, arrays = read_r2d2plus_z(path, names=["state"], verify=verify)
        flat = arrays["state"].reshape(-1)
    else:
        flat = np.frombuffer(raw, dtype="<f8")

    nx = grid["nx"] // mpi["size_x"]
    ny = grid["ny"] // mpi["size_y"]
    nz = grid["nz"] // mpi["size_z"]
    ncomp = len(RESTART_COMPONENTS)
    if flat.size != ncomp * nx * ny * nz:
        raise ValueError(
            f"{path}: {flat.size} values, expected {ncomp}x{nx}x{ny}x{nz} "
            "(state_subdir or the decomposition may be different)"
        )
    qq = flat.reshape((ncomp, nx, ny, nz))
    state = {name: qq[m] for m, name in enumerate(RESTART_COMPONENTS)}
    return state, meta
