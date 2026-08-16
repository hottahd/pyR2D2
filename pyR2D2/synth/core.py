"""スペクトル合成 (3D スナップショットの後処理)。

gray または多群で作った R2D2 のスナップショットから、視線に沿って
輻射輸送を形式解で解き、高分解能スペクトルを作る。

    3D (T, P, v_x) の柱  →  各波長で tau を積む  →  I_lambda  →  水平平均
                                                    →  装置プロファイルで畳む

前提と近似:
    * LTE (S_lambda = B_lambda(T))。散乱は吸収として扱う
    * 不透明度は (T,P) グリッド上に事前計算した高分解能テーブルから補間
      (`scripts/build_opacity_table.py`)
    * 視線速度による Doppler シフトは、等分解能格子上の**添字のずらし**で入れる
    * まず mu = 1 (disk center) の鉛直光線。傾いた光線は将来の拡張
      (箱の中を斜めに補間して辿る必要がある)

**マイクロ乱流を入れないこと**: 3D の速度場が本来の乱流を含んでいるので、
1D モデルで使う vturb を足すと二重計上になる。テーブルは vturb=0 で作る。
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from . import constants as cst
from .means import planck_lambda


@dataclass
class OpacityTable:
    """(T, P) グリッド上の高分解能不透明度テーブル。"""

    T: np.ndarray          # (nT,)
    P: np.ndarray          # (nP,)
    lam: np.ndarray        # (nlam,) [Angstrom], 等分解能 (log 等間隔)
    logk: np.ndarray       # (nT, nP, nlam) log10 kappa [cm^2/g]
    meta: dict = None

    @classmethod
    def load(cls, path):
        import json

        import h5py

        with h5py.File(path, "r") as f:
            return cls(T=f["T"][:], P=f["P"][:], lam=f["lam"][:],
                       logk=f["logk"][:],
                       meta=json.loads(f.attrs.get("meta", "{}")))

    @property
    def resolving_power(self) -> float:
        return 1.0 / (self.lam[1] / self.lam[0] - 1.0)

    def interpolate_column(self, T: np.ndarray, P: np.ndarray) -> np.ndarray:
        """柱の各深さ点での log10 kappa(lambda) を返す (ndepth, nlam)。

        (log T, log P) の双線形補間。範囲外は端で打ち切る。
        """
        lt = np.log10(np.clip(T, self.T[0], self.T[-1]))
        lp = np.log10(np.clip(P, self.P[0], self.P[-1]))
        Tg, Pg = np.log10(self.T), np.log10(self.P)

        i = np.clip(np.searchsorted(Tg, lt), 1, len(Tg) - 1)
        j = np.clip(np.searchsorted(Pg, lp), 1, len(Pg) - 1)
        wt = ((lt - Tg[i - 1]) / (Tg[i] - Tg[i - 1]))[:, None]
        wp = ((lp - Pg[j - 1]) / (Pg[j] - Pg[j - 1]))[:, None]

        return ((1 - wt) * (1 - wp) * self.logk[i - 1, j - 1]
                + wt * (1 - wp) * self.logk[i, j - 1]
                + (1 - wt) * wp * self.logk[i - 1, j]
                + wt * wp * self.logk[i, j])


def doppler_shift_rows(logk: np.ndarray, v_los: np.ndarray,
                       resolving_power: float) -> np.ndarray:
    """各深さの不透明度を視線速度で Doppler シフトする。

    等分解能格子 lambda_n = lambda_0 (1+1/R)^n の上では、波長を (1 - v/c) 倍する
    ことは **添字を一定量ずらす**ことと同じ:

        shift = ln(1 - v/c) / ln(1 + 1/R)   [格子点数]

    観測者の方向 (+x) へ動く要素 (v_los > 0) は青方偏移するので、
    観測波長 lambda で見える不透明度は共動系の lambda(1 - v/c) の値。

    Parameters
    ----------
    logk : (ndepth, nlam)
    v_los : (ndepth,)  視線速度 [cm/s]、観測者向きが正
    """
    nlam = logk.shape[1]
    dln = np.log1p(1.0 / resolving_power)
    shift = np.log1p(-np.asarray(v_los) / cst.c_light) / dln  # 格子点数
    out = np.empty_like(logk)
    idx = np.arange(nlam)
    for k in range(logk.shape[0]):
        out[k] = np.interp(idx + shift[k], idx, logk[k],
                           left=logk[k, 0], right=logk[k, -1])
    return out


def formal_solution(alpha: np.ndarray, source: np.ndarray, dx: np.ndarray,
                    mu: float = 1.0) -> np.ndarray:
    """柱に沿った形式解。上端から出てくる強度を返す。

    深さ配列は **上端が先頭** (alpha[0] が最上層) であること。
    source function が光学的厚みについて区分線形であるとして厳密に積分する。

    Parameters
    ----------
    alpha : (ndepth, nlam)   吸収係数 rho*kappa [cm^-1]
    source : (ndepth, nlam)  source function (LTE なら B_lambda)
    dx : (ndepth-1,)         隣り合う深さ点の間隔 [cm] (正)
    mu : float               cos(視線角)

    Returns
    -------
    I : (nlam,)  上端での強度
    """
    ndepth = alpha.shape[0]
    # 最深部は tau >> 1 とみなして S で初期化 (拡散極限)
    I = source[-1].copy()
    for k in range(ndepth - 2, -1, -1):
        d = 0.5 * (alpha[k] + alpha[k + 1]) * dx[k] / mu
        d = np.maximum(d, 1e-12)
        ex = np.exp(-d)
        e0 = 1.0 - ex
        # (e0 - d*ex)/d は d -> 0 で d/2。小さい d では級数展開で桁落ちを防ぐ
        small = d < 1e-4
        w1 = np.where(small, 0.5 * d, (e0 - d * ex) / d)
        I = I * ex + source[k] * e0 + (source[k + 1] - source[k]) * w1
    return I


def synth_column(T: np.ndarray, P: np.ndarray, rho: np.ndarray,
                 v_los: np.ndarray, x: np.ndarray, table: OpacityTable,
                 mu: float = 1.0) -> np.ndarray:
    """1 本の柱から出射強度スペクトルを作る。

    T, P, rho, v_los, x はいずれも **上端が先頭** で並んでいること。
    x は高さ [cm] (上端が大きい)。
    """
    logk = table.interpolate_column(T, P)
    if np.any(v_los != 0.0):
        logk = doppler_shift_rows(logk, v_los, table.resolving_power)
    alpha = (10.0 ** logk) * rho[:, None]
    S = np.array([planck_lambda(table.lam, float(t)) for t in T])
    dx = np.abs(np.diff(x))
    return formal_solution(alpha, S, dx, mu=mu)


def instrument_profile(lam: np.ndarray, spec: np.ndarray,
                       resolving_power: float) -> np.ndarray:
    """装置プロファイル (Gauss) で畳み込んで分解能を落とす。

    等分解能格子なので、畳み込みは添字空間で一定幅の Gauss になる。
    FWHM = lambda/R に対応する標準偏差を格子点数で表す。
    """
    R_grid = 1.0 / (lam[1] / lam[0] - 1.0)
    fwhm_pix = R_grid / resolving_power
    sigma = fwhm_pix / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    if sigma < 0.3:
        return spec.copy()
    n = int(np.ceil(4 * sigma))
    t = np.arange(-n, n + 1)
    kern = np.exp(-0.5 * (t / sigma) ** 2)
    kern /= kern.sum()
    return np.convolve(spec, kern, mode="same")


def bisector(lam: np.ndarray, spec: np.ndarray, n_level: int = 20,
             lam0: float = None, half_width: float = 0.5):
    """吸収線のバイセクタを求める。

    Parameters
    ----------
    lam0 : float
        線中心の目安 [Angstrom]。省略時はスペクトルの最小値の位置。
    half_width : float
        線の左右に取る幅 [Angstrom]。

    Returns
    -------
    depth : (n_level,)   規格化強度 (0=コア, 1=連続光)
    lam_b : (n_level,)   バイセクタの波長 [Angstrom]
    """
    if lam0 is None:
        lam0 = lam[np.argmin(spec)]
    m = np.abs(lam - lam0) < half_width
    ll, ss = lam[m], spec[m]
    cont = np.percentile(ss, 95)
    core_i = int(np.argmin(ss))
    core = ss[core_i]

    levels = np.linspace(core + 0.05 * (cont - core), cont - 0.05 * (cont - core),
                         n_level)
    lam_b = np.full(n_level, np.nan)
    for i, lv in enumerate(levels):
        # 左側 (青) と右側 (赤) で交差点を線形補間
        try:
            l_blue = np.interp(lv, ss[:core_i + 1][::-1], ll[:core_i + 1][::-1])
            l_red = np.interp(lv, ss[core_i:], ll[core_i:])
            lam_b[i] = 0.5 * (l_blue + l_red)
        except Exception:
            pass
    return (levels - core) / (cont - core), lam_b
