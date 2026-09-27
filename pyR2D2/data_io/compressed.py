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
    if verify is not False and checksums is not None and rank < len(checksums):
        digest = _xxh3_hex(raw)
        if digest is None and verify:
            raise ImportError("verify=True needs the xxhash module")
        if digest is not None and digest != checksums[rank]:
            raise ValueError(f"{path}: XXH3 mismatch ({digest} != {checksums[rank]})")

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
