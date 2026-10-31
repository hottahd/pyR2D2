"""
EUV・軟 X 線の強度を emission measure (EM) と応答関数から作る (R2D2plus DEC-671)。

R2D2plus の量の並び 2 のランは ``remap/vl/vl_dem.dac`` に、視線ごとの EM = Σ n_e² dl [cm^-5] を
log10 T の箱 ([5, 8)、0.05 dex、60 箱) に分けて書く (:py:meth:`pyR2D2.OnTheFly.read_dem`)。
強度は EM に装置の温度応答 R(T) [DN cm^5 s^-1 pix^-1] を掛けて足したもの::

    I = Σ_箱 EM_箱 × (箱の中の R の平均)   [DN s^-1 pix^-1]

使い方::

    d = pyR2D2.Data(run)
    d.vc.read(n)                                   # em_x (鉛直の視線), em_y, em_z
    img = pyR2D2.synth.euv.intensity(d.vc.em_x, "aia")     # {"A94": (jx, kx), ...}
    xrt = pyR2D2.synth.euv.intensity(d.vc.em_x, "xrt", channels=["Al-poly"])

応答関数は :class:`ResponseTable` で、表を差し替えられる (``intensity(em, "my_table.txt")`` や
``ResponseTable`` を直接渡す)。同梱の表 (``pyR2D2/synth/data/``) は SSW の公開ファイルから
:func:`aia_response_from_ssw` と :func:`xrt_response_from_ssw` で作ったもの:

- ``aia``: ``aia_get_response(/temp, /dn, /evenorm, /chiantifix)`` (V10、時間変化の補正なし) と同じ手順。
  94・131・171・193・211・335 Å。304 Å は光学的に厚いので入れない。
- ``xrt``: ``make_xrt_temp_resp(make_xrt_wave_resp(contam_time=...), /chianti_default)`` と同じ手順。
  CCD と焦点面フィルターの汚染 (DEHP) の厚さは表の作った日付のもの (表の頭の ``meta`` に記録)。
  観測と比べるときは、その日付で作り直す (``python -m pyR2D2.synth.euv xrt --date ...``)。

EM が無いラン (古いラン、球座標など) は :func:`em_from_state` で 3D 出力から作り直せる
(EOS 表に ``ne`` が要る。R2D2plus-input の ``stratification_eos/add_electron_density.py``)。
"""

from __future__ import annotations

import argparse
import datetime as _dt
import json
import sys
from pathlib import Path

import numpy as np

#: R2D2plus の vl_dem.dac の箱の境界 (log10 T、[5, 8) を 0.05 dex、60 箱)
DEM_LOGTE_EDGES = 5.0 + 0.05 * np.arange(61)

#: AIA の光学的に薄いチャンネル (304 Å は除く)
AIA_CHANNELS = ("A94", "A131", "A171", "A193", "A211", "A335")

#: XRT の単独フィルター (SSW の xrt_channels の名前)。組み合わせ ("Al-poly/Ti-poly" など) も作れる
XRT_CHANNELS = ("Al-mesh", "Al-poly", "C-poly", "Ti-poly", "Be-thin", "Be-med", "Al-med", "Al-thick", "Be-thick")

_DATA_DIR = Path(__file__).resolve().parent / "data"
#: 同梱の表 (名前 → ファイル)
BUNDLED = {"aia": "aia_response.txt", "xrt": "xrt_response.txt"}

_FORMAT_LINE = "# pyR2D2-response-table 1"


class ResponseTable:
    """
    温度応答の表 R(log10 T) をチャンネルごとに持つ。

    Parameters
    ----------
    logt : array (n,)
        log10 T [K] (増加順)
    responses : dict of str to array (n,)
        チャンネル名 → 応答 (単位は ``units``)
    units : str
        応答の単位。既定は ``"DN cm^5 s^-1 pix^-1"``
    meta : dict, optional
        出どころなどの記録 (保存すると表の頭に JSON で残る)

    表の形式 (テキスト)::

        # pyR2D2-response-table 1
        # units: DN cm^5 s^-1 pix^-1
        # meta: {"source": ..., ...}
        # logT A94 A131 ...
        5.00 1.2e-28 ...

    ``#`` 以外の行は空白区切りの数。``.npz`` (キー ``logt`` と各チャンネル、任意で ``units``・``meta``) も読める。
    """

    def __init__(self, logt, responses, units="DN cm^5 s^-1 pix^-1", meta=None):
        self.logt = np.asarray(logt, dtype=np.float64)
        if self.logt.ndim != 1 or np.any(np.diff(self.logt) <= 0):
            raise ValueError("logt must be one-dimensional and increasing")
        self.responses = {}
        for ch, r in responses.items():
            r = np.asarray(r, dtype=np.float64)
            if r.shape != self.logt.shape:
                raise ValueError(f"{ch}: {r.shape} values for {self.logt.size} temperatures")
            self.responses[str(ch)] = r
        self.units = units
        self.meta = dict(meta or {})

    @property
    def channels(self):
        """チャンネル名の list"""
        return list(self.responses)

    def __getitem__(self, channel):
        return self.responses[channel]

    def __repr__(self):
        return (
            f"ResponseTable({self.channels}, logT {self.logt[0]:.2f}..{self.logt[-1]:.2f}, "
            f"{self.units!r}, source={self.meta.get('source', '?')!r})"
        )

    def bin_mean(self, channel, logte_edges=DEM_LOGTE_EDGES, nsub=20):
        """
        log10 T の箱ごとの応答の平均 (箱の中を log10 T で等分した ``nsub`` 点の平均、表の外は 0)。

        箱の中の DEM の形は分からないので、log10 T について一様と見なしている。
        """
        edges = np.asarray(logte_edges, dtype=np.float64)
        frac = (np.arange(nsub) + 0.5) / nsub
        pts = edges[:-1, None] + np.diff(edges)[:, None] * frac[None, :]
        r = np.interp(pts, self.logt, self.responses[channel], left=0.0, right=0.0)
        return r.mean(axis=1)

    def save(self, path):
        """テキスト (``.npz`` で終わる名前なら npz) で保存する。"""
        path = Path(path)
        if path.suffix == ".npz":
            np.savez(path, logt=self.logt, units=self.units, meta=json.dumps(self.meta), **self.responses)
            return
        for ch in self.channels:
            if any(c.isspace() for c in ch):
                raise ValueError(f"channel name {ch!r} contains a space")
        lines = [
            _FORMAT_LINE,
            f"# units: {self.units}",
            f"# meta: {json.dumps(self.meta, ensure_ascii=False)}",
            "# logT " + " ".join(self.channels),
        ]
        for i, t in enumerate(self.logt):
            lines.append(f"{t:.4f} " + " ".join(f"{self.responses[ch][i]:.6e}" for ch in self.channels))
        path.write_text("\n".join(lines) + "\n")

    @classmethod
    def load(cls, path):
        """テキストか ``.npz`` の表を読む。"""
        path = Path(path)
        if path.suffix == ".npz":
            with np.load(path, allow_pickle=False) as z:
                units = str(z["units"]) if "units" in z.files else "DN cm^5 s^-1 pix^-1"
                meta = json.loads(str(z["meta"])) if "meta" in z.files else {}
                responses = {k: z[k] for k in z.files if k not in ("logt", "units", "meta")}
                return cls(z["logt"], responses, units, meta)
        units, meta, names, rows = "DN cm^5 s^-1 pix^-1", {}, None, []
        for line in path.read_text().splitlines():
            s = line.strip()
            if not s:
                continue
            if s.startswith("#"):
                body = s[1:].strip()
                if body.startswith("units:"):
                    units = body[len("units:"):].strip()
                elif body.startswith("meta:"):
                    meta = json.loads(body[len("meta:"):].strip())
                elif body.split()[:1] == ["logT"]:
                    names = body.split()[1:]
                continue
            rows.append([float(v) for v in s.split()])
        if names is None:
            raise ValueError(f"{path}: no '# logT <channels>' header line")
        a = np.array(rows, dtype=np.float64)
        if a.ndim != 2 or a.shape[1] != len(names) + 1:
            raise ValueError(f"{path}: rows do not have {len(names) + 1} columns")
        return cls(a[:, 0], {ch: a[:, i + 1] for i, ch in enumerate(names)}, units, meta)


def load_response(response="aia"):
    """
    応答の表を返す。``"aia"``・``"xrt"`` は同梱の表、文字列や Path はファイル、
    :class:`ResponseTable` はそのまま。
    """
    if isinstance(response, ResponseTable):
        return response
    key = str(response)
    if key in BUNDLED:
        return ResponseTable.load(_DATA_DIR / BUNDLED[key])
    return ResponseTable.load(response)


def intensity(em, response="aia", channels=None, logte_edges=DEM_LOGTE_EDGES, nsub=20):
    """
    箱ごとの EM から強度を作る。

    Parameters
    ----------
    em : array (..., nbin)
        log10 T の箱ごとの EM [cm^-5] (最後の軸が箱)。``d.vc.em_x`` など
    response : str, Path or ResponseTable
        応答の表 (:func:`load_response`)。既定は同梱の AIA
    channels : list of str, optional
        作るチャンネル。既定は表の全部
    logte_edges : array (nbin + 1,)
        箱の境界。既定は R2D2plus の vl_dem.dac のもの

    Returns
    -------
    dict of str to array (...)
        チャンネル名 → 強度 (単位は応答 × cm^-5。既定の表では DN s^-1 pix^-1)
    """
    table = load_response(response)
    em = np.asarray(em)
    edges = np.asarray(logte_edges, dtype=np.float64)
    if em.shape[-1] != edges.size - 1:
        raise ValueError(f"em has {em.shape[-1]} bins on its last axis, logte_edges gives {edges.size - 1}")
    if channels is None:
        channels = table.channels
    out = {}
    for ch in channels:
        r = table.bin_mean(ch, edges, nsub)
        out[ch] = np.tensordot(em.astype(np.float64), r, axes=([-1], [0]))
    return out


def em_from_state(data, n, los="x", logte_edges=DEM_LOGTE_EDGES):
    """
    3D 出力 (``remap/qq``) から箱ごとの EM を作り直す (R2D2plus の vl_dem.dac と同じ定義)。

    n_e と T は EOS 表 (ρ, s) から引く (表と線形の合成はしない。R2D2plus と同じ)。
    表に ``ne`` が要る (``data.p.read_eos_table`` が読んだ表に ``ne`` があると ``log_ne_e`` ができる)。
    箱の外の温度は数えない。セルの幅は margin 込みの座標 (``p.xg``・``yg``・``zg``) の
    中心差分 0.5 (c[i+1] - c[i-1]) で、物理セルの両端も同じ (R2D2plus の write_hist_dem と同じ)。

    Parameters
    ----------
    data : pyR2D2.Data
    n : int
        3D 出力の番号 (``nd``。vl の番号 ``nd_vl`` とは違うことがある)
    los : {"x", "y", "z"}
        視線の向き (pyR2D2 の軸、x が鉛直)。戻り値の形は x: (jx, kx, nbin)、y: (ix, kx, nbin)、
        z: (ix, jx, nbin) で、:py:meth:`pyR2D2.OnTheFly.read_dem` の ``em_x``・``em_y``・``em_z`` と同じ
    logte_edges : array
        箱の境界 (等間隔でなくてよい)

    Returns
    -------
    numpy.ndarray (float64)
    """
    from pyR2D2.util.util import eos_table

    p = data.p
    if "log_ne_e" not in p.__dict__:
        raise ValueError(
            "the EoS table has no electron density ne; read a table with ne via "
            "data.p.read_eos_table(path) (R2D2plus-input stratification_eos/add_electron_density.py)"
        )
    if p.geometry != "Cartesian":
        raise ValueError(f"em_from_state supports Cartesian runs only (this run is {p.geometry})")
    data.qf.read(n, keys=["ro", "se"])
    ro = data.qf.ro.astype(np.float64) + np.asarray(p.ro0, dtype=np.float64)[:, None, None]
    se = data.qf.se.astype(np.float64) + np.asarray(p.se0, dtype=np.float64)[:, None, None]
    te = eos_table(data, ro, se, "te")
    ne = eos_table(data, ro, se, "ne")
    del ro, se

    edges = np.asarray(logte_edges, dtype=np.float64)
    nb = edges.size - 1
    b = np.searchsorted(edges, np.log10(te), side="right") - 1
    ok = (b >= 0) & (b < nb)
    ix, jx, kx = te.shape
    del te
    ne2 = ne * ne
    del ne

    if los == "x":
        dl = _cell_width(p.xg, p.margin, ix)[:, None, None]
        cell = (np.arange(jx)[:, None] * kx + np.arange(kx)[None, :])[None, :, :]
        shape = (jx, kx)
    elif los == "y":
        dl = _cell_width(p.yg, p.margin, jx)[None, :, None]
        cell = (np.arange(ix)[:, None] * kx + np.arange(kx)[None, :])[:, None, :]
        shape = (ix, kx)
    elif los == "z":
        dl = _cell_width(p.zg, p.margin, kx)[None, None, :]
        cell = (np.arange(ix)[:, None] * jx + np.arange(jx)[None, :])[:, :, None]
        shape = (ix, jx)
    else:
        raise ValueError(f"los must be 'x', 'y' or 'z', not {los!r}")
    w = ne2 * dl
    idx = np.broadcast_to(cell, b.shape) * nb + b
    em = np.bincount(idx[ok], weights=w[ok], minlength=shape[0] * shape[1] * nb)
    return em.reshape(shape + (nb,))


def _cell_width(coord_with_margin, margin, n):
    """
    セルの幅 0.5 (c[i+1] - c[i-1]) を margin 込みの座標から作る (R2D2plus の write_hist_dem と同じ)。
    物理セルの両端も margin の座標を使った中心差分になる。
    """
    c = np.asarray(coord_with_margin, dtype=np.float64)
    if c.size != n + 2 * margin or margin < 1:
        raise ValueError(f"coordinate with margin has {c.size} points, expected {n} + 2 x {margin}")
    i = np.arange(n) + margin
    return 0.5 * (c[i + 1] - c[i - 1])


# ---------------------------------------------------------------------------
# 応答の表を SSW の公開ファイルから作る
# ---------------------------------------------------------------------------

#: SSW の AIA の応答のファイルの置き場
SSW_AIA_URL = "https://hesperia.gsfc.nasa.gov/ssw/sdo/aia/response/"
#: SSW の XRT の応答のファイルの置き場
SSW_XRT_URL = "https://hesperia.gsfc.nasa.gov/ssw/hinode/xrt/idl/response/"


def _read_genx(path):
    from sunpy.io.special import read_genx

    return read_genx(str(path))


def _aia_response_table_rows(path):
    """aia_V*_response_table.txt を行の dict の list で読む。"""
    lines = [s for s in Path(path).read_text().splitlines() if s.strip()]
    head = lines[0].split()
    rows = []
    for s in lines[1:]:
        w = s.split()
        rows.append(dict(zip(head, w)))
    return rows


def _parse_utc(s):
    s = s.rstrip("Z")
    for fmt in ("%Y-%m-%dT%H:%M:%S.%f", "%Y-%m-%dT%H:%M:%S", "%Y-%m-%d"):
        try:
            return _dt.datetime.strptime(s, fmt)
        except ValueError:
            pass
    raise ValueError(f"cannot read the time {s!r}")


def aia_response_from_ssw(
    ssw_dir,
    channels=AIA_CHANNELS,
    evenorm=True,
    chiantifix=True,
    date=None,
    version="V10",
    response_table=None,
    blend=True,
):
    """
    AIA の温度応答を SSW の公開ファイルから作る (``aia_get_response(/temp, /dn, ...)`` と同じ手順)。

    要るファイル (``ssw_dir`` に置く、:py:data:`SSW_AIA_URL` から取れる)::

        aia_V9_all_fullinst.genx      有効面積 (V10 も V8 と同じ中身。IDL の V10 も V9 の名前を読む)
        aia_V9_fullemiss.genx         CHIANTI の放射率 (790 MB)
        aia_V9_chiantifix.genx        94・131 Å の欠けた線の経験的な補正 (chiantifix のとき)
        aia_V10_20201119_190000_response_table.txt   EVE の規格化と時間変化 (evenorm・date のとき)

    手順 (SSW の aia_bp_blend_channels・aia_bp_parse_effarea・aia_bp_corrections・aia_bp_make_tresp・
    aia_bp_parse_tresp):

    1. ``blend``: 94←304、131←335、335←131 のクロストークを有効面積に足す (IDL の ``/noblend`` は ``blend=False``)
    2. DN へ: × max(12398/λ × elecperev, 1) / elecperdn
    3. ``evenorm``: × (表の最初の期間の有効面積) / (有効波長での有効面積)
    4. ``date``: × その日付の期間の比 × (1 + p1 dt + p2 dt² + p3 dt³) (dt は期間の始めからの日数)
    5. R(T) = Σ_λ 有効面積(λ) × 放射率(λ, T) × platescale × Δλ
    6. ``chiantifix``: 94・131 Å に経験的な差 (``date`` のときは 4 の比を掛けたもの) を足す

    Returns
    -------
    ResponseTable
    """
    ssw_dir = Path(ssw_dir)
    inst = _read_genx(ssw_dir / "aia_V9_all_fullinst.genx")

    def full(ch):
        return inst[f"{ch}_FULL"]

    # 1. クロストーク (aia_bp_blend_channels)。元の (混ぜる前の) 成分から作る
    cross_of = {"A94": "A304", "A131": "A335", "A304": "A94", "A335": "A131"}
    ea_cm2 = {}
    for ch in set(channels) | {"A304"}:
        f = full(ch)
        ea = np.asarray(f["EFFAREA"], dtype=np.float64)
        if blend and ch in cross_of:
            c = full(cross_of[ch])
            ea = ea + (
                np.asarray(c["ENT_FILTER"], np.float64) * c["PRIMARY"] * c["SECONDARY"] * c["CONTAM"]
                * c["CCD"] * c["GEOAREA"] * np.asarray(f["FP_FILTER"], np.float64)
            )
        ea_cm2[ch] = ea
    wave = np.asarray(inst[channels[0]]["WAVE"], dtype=np.float64)
    platescale = float(inst[channels[0]]["PLATESCALE"])

    # 3, 4. EVE の規格化と時間変化 (aia_bp_corrections)
    scale = {ch: 1.0 for ch in channels}
    tscale = {ch: 1.0 for ch in channels}
    table_name = None
    if evenorm or date is not None:
        if response_table is None:
            cands = sorted(ssw_dir.glob(f"aia_{version}_*_response_table.txt"))
            if not cands:
                raise FileNotFoundError(f"no aia_{version}_*_response_table.txt in {ssw_dir}")
            response_table = cands[-1]
        table_name = Path(response_table).name
        rows = _aia_response_table_rows(response_table)
        when = None if date is None else _parse_utc(str(date))
        for ch in channels:
            mine = [r for r in rows if r["WAVELNTH"] == ch[1:] and r["WAVE_STR"].endswith("_THIN")]
            mine.sort(key=lambda r: _parse_utc(r["T_START"]))
            first = mine[0]
            if evenorm:
                nom = np.interp(float(first["EFF_WVLN"]), wave, ea_cm2[ch])
                scale[ch] = float(first["EFF_AREA"]) / nom
            if when is not None:
                pre = [r for r in mine if _parse_utc(r["T_START"]) <= when]
                if not pre:
                    raise ValueError(f"{date} is before the first epoch of {ch}")
                use = pre[-1]
                dt = (when - _parse_utc(use["T_START"])).total_seconds() / 86400.0
                p1, p2, p3 = (float(use[k]) for k in ("EFFA_P1", "EFFA_P2", "EFFA_P3"))
                tscale[ch] = float(use["EFF_AREA"]) / float(first["EFF_AREA"]) * (1 + p1 * dt + p2 * dt**2 + p3 * dt**3)

    # 5. 放射率との積分 (aia_bp_make_tresp)
    em = _read_genx(ssw_dir / "aia_V9_fullemiss.genx")["TOTAL"]
    ewave = np.asarray(em["WAVE"], dtype=np.float64)
    logt = np.asarray(em["LOGTE"], dtype=np.float64)
    emiss = np.asarray(em["EMISSIVITY"], dtype=np.float64)
    if emiss.shape == (logt.size, ewave.size):
        emiss = emiss.T
    dw = ewave[1] - ewave[0]
    responses = {}
    for ch in channels:
        f = full(ch)
        evperphot = 12398.0 / wave
        elecperphot = np.maximum(evperphot * float(f["ELECPEREV"]), 1.0)
        ea_dn = ea_cm2[ch] * elecperphot / float(f["ELECPERDN"]) * scale[ch] * tscale[ch]
        iresp = np.interp(ewave, wave, ea_dn)
        iresp[(ewave < wave.min()) | (ewave > wave.max())] = 0.0
        responses[ch] = (np.maximum(iresp, 0.0) @ emiss) * platescale * dw

    # 6. chiantifix (aia_bp_parse_tresp)
    if chiantifix:
        fix = _read_genx(ssw_dir / "aia_V9_chiantifix.genx")
        fix_ch = [str(c).strip().upper() for c in np.atleast_1d(fix["CHANNELS"])]
        diff = np.asarray(fix["EMPIRICAL_MINUS_RAW"], dtype=np.float64)
        flogt = np.asarray(fix["LOGTE"], dtype=np.float64)
        if diff.shape[0] != flogt.size:
            diff = diff.T
        if flogt.size != logt.size or np.max(np.abs(flogt - logt)) > 1e-6:
            raise ValueError("chiantifix logT grid differs from the emissivity grid")
        for ch in channels:
            if ch.upper() in fix_ch:
                responses[ch] = responses[ch] + diff[:, fix_ch.index(ch.upper())] * tscale[ch]

    meta = {
        "source": "SSW AIA response (aia_get_response /temp /dn"
        + (" /evenorm" if evenorm else "")
        + (" /chiantifix" if chiantifix else "")
        + (f" timedepend_date={date}" if date is not None else "")
        + ("" if blend else " /noblend")
        + f"), version {version}",
        "files": ["aia_V9_all_fullinst.genx", "aia_V9_fullemiss.genx"]
        + (["aia_V9_chiantifix.genx"] if chiantifix else [])
        + ([table_name] if table_name else []),
        "url": SSW_AIA_URL,
        "made_by": "pyR2D2.synth.euv.aia_response_from_ssw",
        "made_on": _dt.date.today().isoformat(),
        "evenorm": bool(evenorm),
        "chiantifix": bool(chiantifix),
        "blend": bool(blend),
        "date": None if date is None else str(date),
        "eve_scale": {ch: scale[ch] for ch in channels},
        "time_scale": {ch: tscale[ch] for ch in channels},
    }
    return ResponseTable(logt, responses, "DN cm^5 s^-1 pix^-1", meta)


def _readsav(path):
    from scipy.io import readsav

    return readsav(str(path))


def _s(v):
    return v.decode() if isinstance(v, bytes) else str(v)


#: SSW xrt_chan2fw の (フィルターホイール 1、2) の位置
_XRT_FW = {
    "Open": (0, 0), "Al-mesh": (0, 1), "Al-poly": (1, 0), "C-poly": (2, 0), "Ti-poly": (0, 2),
    "Be-thin": (3, 0), "Be-med": (4, 0), "Al-med": (5, 0), "Al-thick": (0, 4), "Be-thick": (0, 5),
    "Al-poly/Al-mesh": (1, 1), "Al-poly/Ti-poly": (1, 2), "Al-poly/Al-thick": (1, 4),
    "Al-poly/Be-thick": (1, 5), "C-poly/Ti-poly": (2, 2), "C-poly/Al-thick": (2, 4),
}

_ANYTIM_EPOCH = _dt.datetime(1979, 1, 1)


def _anytim(when):
    """SSW anytim の秒 (1979-01-01 からの秒、閏秒なし)。"""
    return (when - _ANYTIM_EPOCH).total_seconds()


def xrt_contamination(ssw_dir, channel, date):
    """
    XRT の汚染 (DEHP) の厚さ [Å] (CCD、焦点面フィルター) を SSW xrt_time2contam と同じく作る。

    ``contam/xrt_contam_on_ccd.geny`` と ``contam/xrt_contam_on_filter.geny`` を時刻で線形に補間し、
    四捨五入する。データの最後の時刻より後は作れない (ValueError)。
    """
    ssw_dir = Path(ssw_dir)
    when = _parse_utc(str(date)) if not isinstance(date, _dt.datetime) else date
    t = _anytim(when)
    ccd = _readsav(_find(ssw_dir, "xrt_contam_on_ccd.geny"))
    filt = _readsav(_find(ssw_dir, "xrt_contam_on_filter.geny"))
    tc, thc = np.asarray(ccd["p1"], np.float64), np.asarray(ccd["p2"], np.float64)
    if t > tc.max():
        raise ValueError(
            f"{when} is after the last CCD contamination time in xrt_contam_on_ccd.geny "
            f"({_ANYTIM_EPOCH + _dt.timedelta(seconds=float(tc.max()))}); fetch a newer file"
        )
    c = float(np.floor(np.interp(t, tc, thc) + 0.5))
    tf, thf = np.asarray(filt["p1"], np.float64), np.asarray(filt["p2"], np.float64)  # (6, 2, ntime)
    fw1, fw2 = _XRT_FW[channel]
    f1 = np.interp(t, tf, thf[fw1, 0])
    f2 = np.interp(t, tf, thf[fw2, 1])
    return c, float(np.floor(f1 + 0.5) + np.floor(f2 + 0.5))


def _find(ssw_dir, name):
    for p in (ssw_dir / name, ssw_dir / "contam" / name, ssw_dir / "channels" / name, ssw_dir / "emiss" / name):
        if p.exists():
            return p
    raise FileNotFoundError(f"{name} not found under {ssw_dir}")


def _dehp_transmission(ssw_dir, thick):
    """DEHP の 1 層の透過率 (SSW multifilt と同じ薄膜の式、垂直入射)。波長 [Å] と透過率を返す。"""
    path = _find(ssw_dir, "n_DEHP.txt")
    a = np.loadtxt(path, skiprows=2)
    wave = a[:, 0] * 10.0
    keep = (wave >= 1.0) & (wave <= 400.0)
    wave, delta, beta = wave[keep], a[keep, 1], a[keep, 2]
    if thick == 0:
        return wave, np.ones_like(wave)
    index = (1.0 - delta) + 1j * beta
    kl = 2.0 * np.pi * index / wave * thick
    kl = kl.real + 1j * np.minimum(kl.imag, 300.0)
    m00 = np.cos(kl)
    m01 = -1j * np.sin(kl) / index
    m10 = -1j * index * np.sin(kl)
    m11 = np.cos(kl)
    t = 2.0 / (m00 + m01 + m10 + m11)
    return wave, np.abs(t) ** 2


def xrt_response_from_ssw(
    ssw_dir,
    channels=XRT_CHANNELS,
    date=None,
    contam_thick=None,
    emiss_file="XRT_emiss_model.default_CHIANTI_20240306.geny",
    channels_file="xrt_channels_v0017.geny",
    port="R",
):
    """
    XRT の温度応答を SSW の公開ファイルから作る
    (``make_xrt_temp_resp(make_xrt_wave_resp(contam_time=date), /chianti_default)`` と同じ手順)。

    要るファイル (``ssw_dir`` か、その ``channels/``・``emiss/``・``contam/`` に置く。
    :py:data:`SSW_XRT_URL` から取れる)::

        channels/xrt_channels_v0017.geny         各チャンネルの透過率・鏡・CCD
        emiss/XRT_emiss_model.default_CHIANTI_20240306.geny   CHIANTI 10 (コロナの組成) のスペクトル
        contam/xrt_contam_on_ccd.geny, contam/xrt_contam_on_filter.geny, contam/n_DEHP.txt   汚染

    手順 (SSW の make_xrt_ea・make_xrt_sr・make_xrt_contam・multifilt・make_xrt_temp_resp):

    1. 汚染の厚さ: ``date`` から (:func:`xrt_contamination`)、または ``contam_thick`` = (CCD, フィルター) [Å]
    2. 有効面積 = 透過率 × DEHP(CCD) × DEHP(フィルター) × 開口面積
    3. 分光応答 = 有効面積 × (hc/λ) × (画素 / 焦点距離)² / 3.65 eV / gain (既定は R 側の読み出し)
    4. R(T) = Σ_λ スペクトル(λ, T) × 分光応答(λ) × Δλ

    IDL は分光応答をスペクトルの波長へ 2 次で補間する (``interpol, /quad``) が、ここは線形。

    Returns
    -------
    ResponseTable
    """
    ssw_dir = Path(ssw_dir)
    if (date is None) == (contam_thick is None):
        raise ValueError("give either date or contam_thick")
    chans = _readsav(_find(ssw_dir, channels_file))["p0"]
    names = [_s(v) for v in chans["NAME"]]
    em = _readsav(_find(ssw_dir, emiss_file))["p0"][0]
    swave = np.asarray(em["WAVE"], np.float64)[: int(em["WLENGTH"])]
    temp = np.asarray(em["TEMP"], np.float64)[: int(em["TLENGTH"])]
    spec = np.asarray(em["SPEC"], np.float64)[: int(em["TLENGTH"]), : int(em["WLENGTH"])]  # (T, λ)
    w0 = np.concatenate([[2 * swave[0] - swave[1]], swave, [2 * swave[-1] - swave[-2]]])
    dsw = (w0[2:] - w0[:-2]) / 2.0

    h, c = 4.136e-15, 2.998e18
    responses, thick = {}, {}
    for ch in channels:
        if ch not in names:
            raise ValueError(f"{ch} is not in {channels_file} ({names})")
        rec = chans[names.index(ch)]
        n = int(rec["LENGTH"])
        wave = np.asarray(rec["WAVE"], np.float64)[:n]
        trans = np.asarray(rec["TRANS"], np.float64)[:n]
        if date is not None:
            tc, tf = xrt_contamination(ssw_dir, ch, date)
        else:
            tc, tf = (float(v) for v in contam_thick)
        thick[ch] = (tc, tf)
        cw, ct = _dehp_transmission(ssw_dir, tc)
        _, ft = _dehp_transmission(ssw_dir, tf)
        inside = (wave >= cw.min()) & (wave <= cw.max())
        wave, trans = wave[inside], trans[inside]
        geom, ccd = rec["GEOM"][0], rec["CCD"][0]
        ea = trans * np.interp(wave, cw, ft) * np.interp(wave, cw, ct) * float(geom["APERTURE_AREA"])
        pix = float(ccd["PIXEL_SIZE"]) / 1e4
        gain = float(ccd["GAIN_R"] if port == "R" else ccd["GAIN_L"])
        sr = ea * (h * c / wave) * (pix / float(geom["FOC_LEN"])) ** 2 / float(ccd["EV_PER_EL"]) / gain
        ss = (swave >= wave.min()) & (swave <= wave.max())
        r = np.interp(swave[ss], wave, sr)
        responses[ch] = spec[:, ss] @ (r * dsw[ss])

    meta = {
        "source": "SSW XRT response (make_xrt_wave_resp + make_xrt_temp_resp)",
        "files": [channels_file, emiss_file] + (
            ["xrt_contam_on_ccd.geny", "xrt_contam_on_filter.geny", "n_DEHP.txt"]
        ),
        "url": SSW_XRT_URL,
        "emission_model": _s(em["NAME"]),
        "made_by": "pyR2D2.synth.euv.xrt_response_from_ssw",
        "made_on": _dt.date.today().isoformat(),
        "date": None if date is None else str(date),
        "contam_thick_angstrom": {ch: {"ccd": thick[ch][0], "filter": thick[ch][1]} for ch in channels},
        "readout_port": port,
    }
    return ResponseTable(np.log10(temp), responses, "DN cm^5 s^-1 pix^-1", meta)


def main(argv=None):
    """``python -m pyR2D2.synth.euv aia|xrt --ssw DIR --out FILE ...`` で応答の表を作る。"""
    ap = argparse.ArgumentParser(prog="python -m pyR2D2.synth.euv", description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("aia", help=f"AIA temperature response from SSW files ({SSW_AIA_URL})")
    a.add_argument("--ssw", required=True, help="directory with the SSW AIA response files")
    a.add_argument("--out", required=True, help="output table (.txt or .npz)")
    a.add_argument("--channels", default=",".join(AIA_CHANNELS))
    a.add_argument("--no-evenorm", action="store_true")
    a.add_argument("--no-chiantifix", action="store_true")
    a.add_argument("--no-blend", action="store_true", help="no crosstalk between channels (IDL /noblend)")
    a.add_argument("--date", help="apply the time-dependent degradation at this UTC time (YYYY-MM-DD[THH:MM:SS])")
    x = sub.add_parser("xrt", help=f"XRT temperature response from SSW files ({SSW_XRT_URL})")
    x.add_argument("--ssw", required=True, help="directory with the SSW XRT response files")
    x.add_argument("--out", required=True, help="output table (.txt or .npz)")
    x.add_argument("--channels", default=",".join(XRT_CHANNELS))
    g = x.add_mutually_exclusive_group(required=True)
    g.add_argument("--date", help="contamination at this UTC time (YYYY-MM-DD[THH:MM:SS])")
    g.add_argument("--contam-thick", help="CCD,filter DEHP thickness in Angstrom, e.g. 0,0 for a clean instrument")
    args = ap.parse_args(argv)
    chans = [c for c in args.channels.split(",") if c]
    if args.cmd == "aia":
        table = aia_response_from_ssw(
            args.ssw, chans, evenorm=not args.no_evenorm, chiantifix=not args.no_chiantifix, date=args.date,
            blend=not args.no_blend,
        )
    else:
        thick = None if args.contam_thick is None else [float(v) for v in args.contam_thick.split(",")]
        table = xrt_response_from_ssw(args.ssw, chans, date=args.date, contam_thick=thick)
    table.save(args.out)
    print(f"wrote {args.out}: {table}")
    return 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
