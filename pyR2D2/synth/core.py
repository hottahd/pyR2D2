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

from . import _physics as cst
from ._physics import planck_lambda


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


def _wrap_bilinear(a3, fy, fz):
    """3D 配列 (nx, ny, nz) を各高さで水平方向に双線形補間する。

    `fy[k]`, `fz[k]` は高さ k での**格子単位の**位置 (実数、周期で巻く)。
    返り値は (nx,)。

    **周期境界は添字を剰余で巻く**。`fmod` は負の値に対して負を返すので
    `% n` を使うこと (2026-08-18: 学生さんの C++ 実装で同じ罠を指摘した)。
    """
    nx, ny, nz = a3.shape
    j0 = np.floor(fy).astype(np.int64) % ny
    k0 = np.floor(fz).astype(np.int64) % nz
    j1 = (j0 + 1) % ny
    k1 = (k0 + 1) % nz
    wy = fy - np.floor(fy)
    wz = fz - np.floor(fz)
    i = np.arange(nx)
    return (a3[i, j0, k0] * (1 - wy) * (1 - wz)
            + a3[i, j1, k0] * wy * (1 - wz)
            + a3[i, j0, k1] * (1 - wy) * wz
            + a3[i, j1, k1] * wy * wz)


def _wrap_trilinear(a3, fx, fy, fz):
    """3D 配列を (鉛直, 水平, 水平) の 3 次元線形補間する。

    `fx` は**格子単位の鉛直位置** (実数、範囲内に丸める)、
    `fy`, `fz` は水平位置 (周期で巻く)。返り値は `fx` と同じ長さ。

    **なぜ要るか (2026-08-18)**: 最初は「立方体全体を鉛直に細分してから
    水平補間する」実装にしていたが、これは**光線 1 本ごとに 128x128 柱を
    全部内挿する**ので実用にならなかった (30 分で 1 例も終わらず)。
    光線が通る点だけを補間すれば無駄がない。
    """
    nx, ny, nz = a3.shape
    fx = np.clip(fx, 0.0, nx - 1.0)
    i0 = np.minimum(np.floor(fx).astype(np.int64), nx - 2)
    i1 = i0 + 1
    wx = fx - i0
    j0 = np.floor(fy).astype(np.int64) % ny
    k0 = np.floor(fz).astype(np.int64) % nz
    j1 = (j0 + 1) % ny
    k1 = (k0 + 1) % nz
    wy = fy - np.floor(fy)
    wz = fz - np.floor(fz)

    def plane(i):
        return (a3[i, j0, k0] * (1 - wy) * (1 - wz)
                + a3[i, j1, k0] * wy * (1 - wz)
                + a3[i, j0, k1] * (1 - wy) * wz
                + a3[i, j1, k1] * wy * wz)

    return plane(i0) * (1 - wx) + plane(i1) * wx


def ray_geometry(x, mu, phi, dy, dz, iy0=0.0, iz0=0.0):
    """傾いた光線が各高さで通る水平位置を、**格子単位**で返す。

    Parameters
    ----------
    x : (nx,)   高さ [cm]、**上端が先頭** (降順)
    mu : float  cos(視線角)。1 が視線中心
    phi : float 方位角 [rad]。視線の水平成分の向き
    dy, dz : float  水平格子幅 [cm]
    iy0, iz0 : float  **上端で**光線が通る位置 [格子単位]

    Returns
    -------
    fy, fz : (nx,)  各高さでの水平位置 [格子単位]
    ds : (nx-1,)    隣り合う高さの間の**経路長** [cm]

    **幾何**: 観測者の方向を n = (mu, sqrt(1-mu^2) cos phi,
    sqrt(1-mu^2) sin phi) とする (x が鉛直上向き)。光線は +n 方向へ進んで
    観測者に届くので、上端から深部へ辿るときは -n 方向へ動く:

        r(s) = r_top - s n,   x(s) = x_top - s mu

    したがって高さ x での水平位置は

        y = y_top - (x_top - x) n_y / mu

    経路長は鉛直間隔 Delta x に対し **Delta s = Delta x / mu**。
    """
    x = np.asarray(x, dtype=float)
    st = np.sqrt(max(0.0, 1.0 - mu * mu))
    ny_dir, nz_dir = st * np.cos(phi), st * np.sin(phi)
    drop = (x[0] - x) / mu                      # 経路長 [cm] (上端から)
    fy = iy0 - drop * ny_dir / dy
    fz = iz0 - drop * nz_dir / dz
    ds = np.abs(np.diff(x)) / mu
    return fy, fz, ds


def _ray_sample(T3, P3, rho3, vx3, vy3, vz3, x, dy, dz,
                mu=1.0, phi=0.0, iy0=0.0, iz0=0.0, nsub=1):
    """傾いた光線の標本点を作り、そこでの物理量を返す。

    :func:`synth_ray` と :func:`ray_contribution` の共通部分。
    **視線速度は鉛直成分と水平成分に分けたまま**返す。限界効果 (limb effect)
    が観測の 2 倍強い原因を、`mu vx` と `sqrt(1-mu^2) v_h` のどちらが
    作っているかで切り分けるため (`docs/12` 13.8 節)。

    Returns
    -------
    xs : (m,)       標本点の高さ [cm] (上端が先頭)
    ds : (m-1,)     隣り合う標本点の**真の経路長** [cm]
    T, P, rho : (m,)
    v_vert : (m,)   `mu * vx`                                     [cm/s]
    v_horiz : (m,)  `sqrt(1-mu^2) * (vy cos phi + vz sin phi)`     [cm/s]

    どちらも符号は :func:`doppler_shift_rows` の約束 (観測者向きが正)。
    """
    x = np.asarray(x, dtype=float)
    nx = len(x)

    # 光線上の標本点を**格子単位の鉛直位置**で作る。nsub=1 なら格子点そのもの。
    # **立方体全体を細分してはいけない** (光線 1 本ごとに全柱を内挿すること
    # になり実用にならない。2026-08-18 に一度その実装で詰まった)。
    if nsub > 1:
        fx = np.concatenate([np.arange(nx - 1)[:, None]
                             + np.linspace(0.0, 1.0, nsub, endpoint=False)[None, :]
                             ]).ravel()
        fx = np.append(fx, float(nx - 1))
    else:
        fx = np.arange(nx, dtype=float)
    # 標本点の高さ [cm] (x は等間隔とは限らないので線形内挿で求める)
    xs = np.interp(fx, np.arange(nx), x)

    st = np.sqrt(max(0.0, 1.0 - mu * mu))
    ny_dir, nz_dir = st * np.cos(phi), st * np.sin(phi)
    drop = (x[0] - xs) / mu                     # 上端からの経路長 [cm]
    fy = iy0 - drop * ny_dir / dy
    fz = iz0 - drop * nz_dir / dz
    ds = np.abs(np.diff(xs)) / mu               # 真の経路長 [cm]

    T = np.exp(_wrap_trilinear(np.log(T3), fx, fy, fz))
    P = np.exp(_wrap_trilinear(np.log(P3), fx, fy, fz))
    rho = np.exp(_wrap_trilinear(np.log(rho3), fx, fy, fz))
    vx = _wrap_trilinear(vx3, fx, fy, fz)
    vy = _wrap_trilinear(vy3, fx, fy, fz)
    vz = _wrap_trilinear(vz3, fx, fy, fz)
    v_vert = mu * vx
    v_horiz = st * (vy * np.cos(phi) + vz * np.sin(phi))
    return xs, ds, T, P, rho, v_vert, v_horiz


def synth_ray(T3, P3, rho3, vx3, vy3, vz3, x, dy, dz, table,
              mu=1.0, phi=0.0, iy0=0.0, iz0=0.0, scheme="log", nsub=1,
              wvert=1.0, whoriz=1.0):
    """**傾いた光線**に沿って形式解を解き、出射強度を返す。

    `synth_column` の 3D 版。光線は箱の下端から上端まで**一本の長特性**として
    通し切り、水平の周期性は**強度を側面から入れ直すのではなく座標を巻く**
    ことで扱う。こうすると初期条件の無い面が現れないので**反復が要らない**
    (`docs/12` 6.1 節。R2D2 本体の RTE が反復しているのは、全セルの J が要る +
    MPI で領域分割している + 短特性、という別条件のため)。

    Parameters
    ----------
    T3, P3, rho3, vx3, vy3, vz3 : (nx, ny, nz)
        **上端が先頭** (x[0] が最上層)。vx は鉛直 (上向きが正)。
    x : (nx,)      高さ [cm]、降順
    dy, dz : float 水平格子幅 [cm]
    mu : float     cos(視線角)
    phi : float    方位角 [rad]
    iy0, iz0 : float  上端での光線の水平位置 [格子単位]
    nsub : int
        鉛直 1 格子あたりの分割数。**1 格子進む間に水平へ 1 格子以上ずれると
        構造を飛び越す**ので、その場合は 2 以上にする
        (:func:`substeps_needed` が必要数を返す)。
    wvert, whoriz : float
        視線速度の**鉛直成分・水平成分に掛ける係数** (既定は 1.0 = そのまま)。
        `whoriz=0` にすると水平速度を消した合成ができる。限界効果が
        どちらの成分で決まっているかを切り分ける実験用 (`docs/12` 13.8 節)。
        **物理を変える操作なので、既定値以外を使ったら必ず記録すること。**

    Returns
    -------
    I : (nlam,)  上端での強度

    **視線速度**: v_los = mu vx + sqrt(1-mu^2) (vy cos phi + vz sin phi)。
    観測者向きが正で、`doppler_shift_rows` の約束と一致する。

    **補間**: 正値の量 (T, P, rho) は **log 空間**、速度は線形。
    形式解が「alpha が経路長に対して指数関数」を仮定している以上、
    端点の値も log で補間しないと一貫しない (HI20 eq. A2 と同じ流儀)。
    """
    xs, ds, T, P, rho, v_vert, v_horiz = _ray_sample(
        T3, P3, rho3, vx3, vy3, vz3, x, dy, dz,
        mu=mu, phi=phi, iy0=iy0, iz0=iz0, nsub=nsub)
    v_los = wvert * v_vert + whoriz * v_horiz

    logk = table.interpolate_column(T, P)
    if np.any(v_los != 0.0):
        logk = doppler_shift_rows(logk, v_los, table.resolving_power)
    alpha = (10.0 ** logk) * rho[:, None]
    S = np.array([planck_lambda(table.lam, float(t)) for t in T])
    # ds は既に真の経路長なので、形式解には mu=1 を渡す
    return formal_solution(alpha, S, ds, mu=1.0, scheme=scheme)


def ray_contribution(T3, P3, rho3, vx3, vy3, vz3, x, dy, dz, table,
                     lam_idx, mu=1.0, phi=0.0, iy0=0.0, iz0=0.0, nsub=1,
                     wvert=1.0, whoriz=1.0):
    """傾いた光線に沿った**寄与関数**と、視線速度の 2 成分を返す。

    出射強度は経路長 s について

        I = int CF ds,      CF(s) = S exp(-tau) alpha

    と書ける (`scripts/contribution_function.py` と同じ定義。ただしあちらは
    mu=1 の鉛直柱専用)。**線位置は「CF で重みを付けた視線速度の平均」で
    ほぼ決まる**ので、その重みを使って `mu vx` と `sqrt(1-mu^2) v_h` の
    寄与を分けて出せる (`docs/12` 13.8 節の分解)。

    Parameters
    ----------
    lam_idx : array_like of int
        寄与関数を返す波長の添字 (例: 線コアと連続光の 2 点)。
        **全波長を返すと (深さ x 波長) が大きくなる**ので絞ること。
    その他 : :func:`synth_ray` と同じ

    Returns
    -------
    dict
        ``xs``     (m,)          標本点の高さ [cm]
        ``T``      (m,)          標本点の温度 [K] (下端の S を作るのに要る)
        ``w``      (m,)          CF に掛ける経路長要素 [cm] (sum(CF*w) = I)
        ``tau``    (m, nsel)     上端から積んだ光学的厚み
        ``cf``     (m, nsel)     寄与関数 S exp(-tau) alpha
        ``v_vert`` (m,)          視線速度の鉛直成分 [cm/s]
        ``v_horiz``(m,)          視線速度の水平成分 [cm/s]

    **注意**: tau は :func:`formal_solution` と同じ対数平均で積むが、
    CF の台形積分と log スキームの形式解は厳密には一致しない。
    絶対値ではなく**重みとして**使うこと。
    """
    lam_idx = np.atleast_1d(np.asarray(lam_idx, dtype=int))
    xs, ds, T, P, rho, v_vert, v_horiz = _ray_sample(
        T3, P3, rho3, vx3, vy3, vz3, x, dy, dz,
        mu=mu, phi=phi, iy0=iy0, iz0=iz0, nsub=nsub)
    v_los = wvert * v_vert + whoriz * v_horiz

    logk = table.interpolate_column(T, P)
    if np.any(v_los != 0.0):
        logk = doppler_shift_rows(logk, v_los, table.resolving_power)
    alpha = (10.0 ** logk[:, lam_idx]) * rho[:, None]
    lam_sel = table.lam[lam_idx]
    S = np.array([planck_lambda(lam_sel, float(t)) for t in T])

    # 上端から積んだ光学的厚み (alpha の対数平均。formal_solution と同じ流儀)
    lg = np.log(alpha[1:] / alpha[:-1])
    amean = np.where(np.abs(lg) < 1e-8, 0.5 * (alpha[:-1] + alpha[1:]),
                     (alpha[1:] - alpha[:-1]) / np.where(np.abs(lg) < 1e-8, 1.0, lg))
    tau = np.empty_like(alpha)
    tau[0] = 0.0
    tau[1:] = np.cumsum(amean * ds[:, None], axis=0)

    # 標本点あたりの経路長要素 (両端は半分)
    w = np.empty(len(xs))
    w[0] = 0.5 * ds[0]
    w[-1] = 0.5 * ds[-1]
    w[1:-1] = 0.5 * (ds[:-1] + ds[1:])

    return {"xs": xs, "T": T, "w": w, "tau": tau,
            "cf": S * np.exp(-tau) * alpha,
            "v_vert": v_vert, "v_horiz": v_horiz}


def substeps_needed(x, mu, dy, dz, safety=1.0):
    """鉛直 1 格子あたり何分割すれば構造を飛び越さないかを返す。

    鉛直に Delta x 進む間の水平移動は `Delta x sqrt(1-mu^2)/mu`。
    これが水平格子幅を超えたら分割する。R2D2 本体の RTE が
    `Delta l = min(Delta x/mu_x, Delta y/mu_y, Delta z/mu_z)` で刻むのと
    同じ考え方 (`docs/12` 6.2 節)。
    """
    if mu >= 1.0:
        return 1
    st = np.sqrt(max(0.0, 1.0 - mu * mu))
    dx = np.abs(np.diff(np.asarray(x, dtype=float))).max()
    shift = dx * st / mu
    return max(1, int(np.ceil(safety * shift / min(dy, dz))))


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
