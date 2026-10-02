#!/usr/bin/env python3
"""線コアの強度がどの高さから来ているかを調べる (寄与関数)。

**なぜ要るか**: R2D2 の箱の上端は光球の +700 km にある (Witzke et al. 2024 が
使った MURaM は +1000 km で「温度極小まで覆う」と明記)。上部境界の影響は
下方に数スケールハイト (上層で H_p ~ 120-150 km) 及ぶので、線コアの形成高度が
上端に近すぎると、**線の深さが境界条件で決まってしまう**。
線が観測より深すぎる (0.218 vs 0.342) 原因の候補として外せない。

**やること**: 出射強度 I = int S exp(-tau) dtau の被積分関数

    CF(z) = S(z) exp(-tau(z)) alpha(z)      [幾何高さについての寄与]

を高さの関数として出し、線コアと連続光それぞれについて
「強度の何 % が上端から何 km 以内で作られるか」を測る。

使い方:
    python -m pyR2D2.synth.contribution --table $ODF_PROJECT_ROOT/work/opac_fe6302.h5 \
        --run $ODF_PROJECT_ROOT/run/d002 --n 78
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


from pyR2D2.synth._physics import planck_lambda
from pyR2D2.synth.core import OpacityTable, doppler_shift_rows


def log_mean(a0, a1):
    with np.errstate(divide="ignore", invalid="ignore"):
        lg = np.log(a1 / a0)
    near = ~np.isfinite(lg) | (np.abs(lg) < 1e-8)
    return np.where(near, 0.5 * (a0 + a1), (a1 - a0) / np.where(near, 1.0, lg))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--table", required=True)
    ap.add_argument("--run", default="/scr/a000/c0234hotta/odf/run/d002")
    ap.add_argument("--n", type=int, default=-1)
    ap.add_argument("--stride", type=int, default=16)
    ap.add_argument("--xmin", type=float, default=-0.5)
    ap.add_argument("--lam-core", type=float, default=6303.2436,
                    help="線中心の真空波長 [A] (Fe I 6301.5012 空気)")
    ap.add_argument("--out", default="/scr/a000/c0234hotta/odf/figs/contribution.png")
    args = ap.parse_args()

    import pyR2D2

    table = OpacityTable.load(args.table)
    d = pyR2D2.Data(str(Path(args.run) / "data"))
    p = d.p
    n = int(p.nd) if args.n < 0 else args.n

    x = np.array(p.x)
    hgt = (x - p.rstar) / 1e8
    order = np.argsort(x)[::-1]
    x, hgt = x[order], hgt[order]
    keep = hgt > args.xmin
    x, hgt = x[keep], hgt[keep]

    d.qr.read(n, keys=["te", "pr", "ro", "vx"])
    te = (np.array(d.qr.te) + np.array(p.te0)[:, None, None])[order][keep]
    pr = (np.array(d.qr.pr) + np.array(p.pr0)[:, None, None])[order][keep]
    ro = (np.array(d.qr.ro) + np.array(p.ro0)[:, None, None])[order][keep]
    vx = np.array(d.qr.vx)[order][keep]

    print(f"{Path(args.run).name}: n={n}, 高さ {len(x)} 点 "
          f"({hgt[0]*1e3:+.0f} .. {hgt[-1]*1e3:+.0f} km)")
    print(f"箱の上端は光球の {hgt[0]*1e3:+.0f} km "
          f"(MURaM/Witzke 2024 は +1000 km)")

    # 線中心と連続光の波長添字
    i_core = int(np.argmin(np.abs(table.lam - args.lam_core)))
    # 連続光は線から十分離れたところ (窓の端)
    i_cont = int(np.argmin(np.abs(table.lam - (table.lam[0] + 2.0))))
    print(f"線コア {table.lam[i_core]:.3f} A / 連続光 {table.lam[i_cont]:.3f} A")

    ny, nz = te.shape[1], te.shape[2]
    dx = np.abs(np.diff(x))
    cf_core = np.zeros(len(x))
    cf_cont = np.zeros(len(x))
    ncol = 0
    for iy in range(0, ny, args.stride):
        for iz in range(0, nz, args.stride):
            T, P, R, V = te[:, iy, iz], pr[:, iy, iz], ro[:, iy, iz], vx[:, iy, iz]
            logk = table.interpolate_column(T, P)
            logk = doppler_shift_rows(logk, V, table.resolving_power)
            alpha = (10.0 ** logk) * R[:, None]
            S = np.array([planck_lambda(table.lam, float(t)) for t in T])
            # 上端から積んだ光学的深さ (対数平均)
            dtau = log_mean(alpha[:-1], alpha[1:]) * dx[:, None]
            tau = np.empty_like(alpha)
            tau[0] = 0.0
            tau[1:] = np.cumsum(dtau, axis=0)
            # CF(z) = S exp(-tau) alpha  (幾何高さについて積分すると I になる)
            cf = S * np.exp(-tau) * alpha
            cf_core += cf[:, i_core]
            cf_cont += cf[:, i_cont]
            ncol += 1
    cf_core /= ncol
    cf_cont /= ncol

    # 上端から積算して何 % が来ているか
    w = np.concatenate([[dx[0]], 0.5 * (np.concatenate([dx, [dx[-1]]])[:-1]
                                        + np.concatenate([[dx[0]], dx])[:-1])])
    w = np.gradient(-x)                       # 幾何厚み [cm] (上端が先頭)
    for tag, cf in (("線コア", cf_core), ("連続光", cf_cont)):
        tot = np.sum(cf * w)
        cum = np.cumsum(cf * w) / tot
        h_km = hgt * 1e3
        i50 = int(np.searchsorted(cum, 0.5))
        i90 = int(np.searchsorted(cum, 0.9))
        # 上端から 300 km 以内の寄与
        m300 = h_km > h_km[0] - 300.0
        f300 = np.sum((cf * w)[m300]) / tot
        m150 = h_km > h_km[0] - 150.0
        f150 = np.sum((cf * w)[m150]) / tot
        print()
        print(f"--- {tag} ---")
        print(f"  寄与の中央値の高さ (50%)  {h_km[i50]:+6.0f} km")
        print(f"  90% が来る高さまで        {h_km[i90]:+6.0f} km")
        print(f"  **上端から 150 km 以内の寄与  {f150*100:5.1f} %**")
        print(f"  **上端から 300 km 以内の寄与  {f300*100:5.1f} %**")

    fig, ax = plt.subplots(1, 2, figsize=(11, 4.4))
    h_km = hgt * 1e3
    for a, (cf, tag) in zip(ax, ((cf_core, "line core"), (cf_cont, "continuum"))):
        a.plot(cf / cf.max(), h_km, lw=1.8)
        a.axhline(h_km[0], color="r", ls="--", lw=1.2,
                  label=f"box top {h_km[0]:+.0f} km")
        a.axhline(h_km[0] - 300, color="orange", ls=":", lw=1.2,
                  label="top - 300 km")
        a.set_xlabel("contribution function (normalised)")
        a.set_ylabel("height [km]")
        a.set_title(tag)
        a.legend(fontsize=8)
        a.grid(alpha=0.3)
    fig.tight_layout()
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=140)
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
