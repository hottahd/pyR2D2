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

    **導出** (2026-08-17 に符号を修正。それまで逆だった):
    観測者へ向かう速度 v (> 0) の物質が出す静止波長 lambda_rest の光は
    lambda_obs = lambda_rest (1 - v/c) に見える (青方偏移)。逆に解くと

        lambda_rest = lambda_obs / (1 - v/c) ~ lambda_obs (1 + v/c)

    したがって**観測波長 lambda_obs で見える不透明度は、共動系の
    lambda_obs (1 + v/c) の値**である。等分解能格子
    lambda_n = lambda_0 (1+1/R)^n の上ではこれは添字を一定量ずらすことで、

        shift = ln(1 + v/c) / ln(1 + 1/R)   [格子点数]
        out[n] = logk[n + shift]

    こうすると静止波長 n0 にある線は n = n0 - shift、すなわち v > 0 で
    **短波長側 (青) へ動く**。正しい向きである。

    以前は `ln(1 - v/c)` を使っており、上昇流が赤方偏移していた。
    その結果、対流ブルーシフトが符号ごと反転して合成スペクトルに
    現れていた (実測: v = +2 km/s で +1940 m/s の赤方偏移)。

    Parameters
    ----------
    logk : (ndepth, nlam)
    v_los : (ndepth,)  視線速度 [cm/s]、観測者向きが正
    """
    nlam = logk.shape[1]
    dln = np.log1p(1.0 / resolving_power)
    shift = np.log1p(np.asarray(v_los) / cst.c_light) / dln  # 格子点数
    out = np.empty_like(logk)
    idx = np.arange(nlam)
    for k in range(logk.shape[0]):
        out[k] = np.interp(idx + shift[k], idx, logk[k],
                           left=logk[k, 0], right=logk[k, -1])
    return out


def formal_solution(alpha: np.ndarray, source: np.ndarray, dx: np.ndarray,
                    mu: float = 1.0, scheme: str = "log") -> np.ndarray:
    """柱に沿った形式解。上端から出てくる強度を返す。

    深さ配列は **上端が先頭** (alpha[0] が最上層) であること。

    Parameters
    ----------
    alpha : (ndepth, nlam)   吸収係数 rho*kappa [cm^-1]
    source : (ndepth, nlam)  source function (LTE なら B_lambda)
    dx : (ndepth-1,)         隣り合う深さ点の間隔 [cm] (正)
    mu : float               cos(視線角)
    scheme : {"log", "linear"}
        セル内で alpha と S をどう変化させると仮定するか。既定は "log"。

    **なぜ log か** (Hotta & Iijima 2020; Hotta & Toriumi 2020 Appendix A):
    光球では H^- のせいで kappa が 1 セルの間に桁で変わる。alpha を算術平均
    すると光学的厚みを系統的に過大評価する。ln(alpha) が経路長について線形
    (= alpha は指数関数的に変化) と仮定すると、平均は**対数平均**

        <alpha> = (a1 - a0) / ln(a1/a0)

    になる。同様に S も、tau について線形だと急勾配のところで**負になり得る**
    のに対し、ln(S) が tau について線形なら常に正で、桁で変わる状況に強い。
    このとき source の寄与は解析的に積分できて

        int_0^dt S exp(-(dt - t')) dt' = (S_d - S_u exp(-dt)) / (1 + r/dt),
        r = ln(S_d / S_u)

    となる (u = 上流 = 深い側, d = 下流 = 浅い側)。
    **粗い格子でも精度が落ちにくい**のが利点。

    Returns
    -------
    I : (nlam,)  上端での強度
    """
    if scheme not in ("log", "linear"):
        raise ValueError(f"未知の scheme: {scheme}")
    if scheme == "log":
        # log スキームは alpha, S が正であることを前提にする。
        # 負の密度・圧力 (k-omega フィルタが上層で作ることがある) が混じると
        # log(負) で NaN になる。linear スキームはこれを黙って飲み込んでしまい
        # 気付けないので、ここで明示的に落とす。
        if not (np.all(alpha > 0.0) and np.all(source > 0.0)):
            n_a = int(np.sum(alpha <= 0.0))
            n_s = int(np.sum(source <= 0.0))
            raise ValueError(
                f"log スキームには正の alpha と S が要る "
                f"(alpha <= 0 が {n_a} 点, S <= 0 が {n_s} 点)。"
                f"入力の密度・圧力・温度が負になっていないか確認すること。")
    ndepth = alpha.shape[0]
    # 最深部は tau >> 1 とみなして S で初期化 (拡散極限)
    I = source[-1].copy()
    for k in range(ndepth - 2, -1, -1):
        a_d, a_u = alpha[k], alpha[k + 1]        # d = 浅い側, u = 深い側
        s_d, s_u = source[k], source[k + 1]

        if scheme == "linear":
            d = 0.5 * (a_d + a_u) * dx[k] / mu
        else:
            # alpha の対数平均。a_d ~ a_u では 0/0 になるので算術平均に落とす
            lg = np.log(a_u / a_d)
            d = np.where(np.abs(lg) < 1e-8,
                         0.5 * (a_d + a_u),
                         (a_u - a_d) / np.where(np.abs(lg) < 1e-8, 1.0, lg)
                         ) * dx[k] / mu
        d = np.maximum(d, 1e-12)
        ex = np.exp(-d)

        if scheme == "linear":
            e0 = 1.0 - ex
            # (e0 - d*ex)/d は d -> 0 で d/2。小さい d では級数展開で桁落ちを防ぐ
            w1 = np.where(d < 1e-4, 0.5 * d, (e0 - d * ex) / d)
            I = I * ex + s_d * e0 + (s_u - s_d) * w1
        else:
            # ln(S) が tau について線形。分母 (d + r) が 0 に近いところは
            # 極限値 S_d * d を使う (b + 1 = 0 の場合に相当)
            r = np.log(s_d / s_u)
            den = d + r
            safe = np.abs(den) > 1e-8
            integ = np.where(safe,
                             d * (s_d - s_u * ex) / np.where(safe, den, 1.0),
                             s_d * d)
            I = I * ex + integ
    return I


def synth_column(T: np.ndarray, P: np.ndarray, rho: np.ndarray,
                 v_los: np.ndarray, x: np.ndarray, table: OpacityTable,
                 mu: float = 1.0, scheme: str = "log") -> np.ndarray:
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
    return formal_solution(alpha, S, dx, mu=mu, scheme=scheme)


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
