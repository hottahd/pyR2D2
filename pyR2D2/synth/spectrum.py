#!/usr/bin/env python3
"""R2D2 のスナップショットからスペクトルを合成する (Phase 5)。

    3D (T, P, v_x) → 柱ごとに形式解 → I_lambda(y,z)
        → 水平平均 → 装置プロファイル (R=3e5) で畳む → バイセクタ

使い方:
    # 1) 不透明度テーブルを作る (1 回だけ)
    python scripts/build_opacity_table.py --lam-min 6280 --lam-max 6320 \
        --resolving-power 1e6 --out $ODF_PROJECT_ROOT/work/opac_fe6302.h5

    # 2) 合成
    python scripts/synth_spectrum.py --table $ODF_PROJECT_ROOT/work/opac_fe6302.h5 \
        --run $ODF_PROJECT_ROOT/run/d001 --n -1 --stride 4 --R-inst 3e5

注意:
    * 3D の速度場が乱流を含むので、テーブルは vturb=0 で作ること
    * gray の計算では光球上層の温度が高めに出るため、強い線のコアは系統的にずれる
      (弱〜中程度の線・バイセクタは影響が小さい)。docs/05_roadmap.md 参照
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from odfgen.synth import (OpacityTable, bisector, instrument_profile,
                          synth_column)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--table", required=True, help="高分解能不透明度テーブル (h5)")
    ap.add_argument("--run", default="/scr/a000/c0234hotta/odf/run/d001")
    ap.add_argument("--n", type=int, default=-1, help="出力番号 (-1 で最新)")
    ap.add_argument("--stride", type=int, default=4,
                    help="柱を何点おきに使うか (4 なら 64x64 本)")
    ap.add_argument("--xmin", type=float, default=-1.0,
                    help="使う深さの下端 [Mm] (これより深いところは切る)")
    ap.add_argument("--R-inst", type=float, default=3.0e5,
                    help="装置分解能 (岡山の分光器を想定)")
    ap.add_argument("--out", default="notebooks/r2d2_spectrum.png")
    ap.add_argument("--save-npz", default=None)
    args = ap.parse_args()

    import pyR2D2

    table = OpacityTable.load(args.table)
    print(f"不透明度テーブル: {table.lam[0]:.2f}-{table.lam[-1]:.2f} A, "
          f"{len(table.lam)} 点, R={table.resolving_power:.3g}, "
          f"vturb={table.meta.get('vturb_kms')} km/s")

    d = pyR2D2.Data(str(Path(args.run) / "data"))
    p = d.p
    n = int(p.nd) if args.n < 0 else args.n
    print(f"スナップショット: n={n} (t={n*p.dtout/60:.0f} min)")

    d.qr.read(n, keys=["te", "pr", "ro", "vx"])
    te = np.array(d.qr.te) + np.array(p.te0)[:, None, None]
    pr = np.array(d.qr.pr) + np.array(p.pr0)[:, None, None]
    ro = np.array(d.qr.ro) + np.array(p.ro0)[:, None, None]
    vx = np.array(d.qr.vx)
    x = np.array(p.x)

    # 上端が先頭になるよう反転し、深いところは切る (tau >> 1 で寄与しない)
    order = np.argsort(x)[::-1]
    x, te, pr, ro, vx = x[order], te[order], pr[order], ro[order], vx[order]
    keep = (x - p.rstar) / 1e8 > args.xmin
    x, te, pr, ro, vx = x[keep], te[keep], pr[keep], ro[keep], vx[keep]
    print(f"使う深さ: {len(x)} 点 "
          f"({(x[0]-p.rstar)/1e8:+.2f} .. {(x[-1]-p.rstar)/1e8:+.2f} Mm), "
          f"T = {te[-1].mean():.0f} K (最深部)")

    ny, nz = te.shape[1], te.shape[2]
    ys = range(0, ny, args.stride)
    zs = range(0, nz, args.stride)
    ncol = len(list(ys)) * len(list(zs))
    print(f"柱: {ncol} 本 (stride={args.stride})")

    t0 = time.time()
    acc = np.zeros(len(table.lam))
    imgs = []
    for a, iy in enumerate(range(0, ny, args.stride)):
        row = []
        for iz in range(0, nz, args.stride):
            I = synth_column(te[:, iy, iz], pr[:, iy, iz], ro[:, iy, iz],
                             vx[:, iy, iz], x, table)
            acc += I
            row.append(I)
        imgs.append(row)
        if a == 0:
            dt = time.time() - t0
            print(f"  1 行 ({len(row)} 柱) に {dt:.1f} s → "
                  f"全体 {dt*len(list(range(0,ny,args.stride)))/60:.1f} 分の見込み")
    mean_spec = acc / ncol
    print(f"合成 {time.time()-t0:.1f} s")

    # 連続光で規格化
    cont = np.percentile(mean_spec, 99)
    norm = mean_spec / cont
    degraded = instrument_profile(table.lam, norm, args.R_inst)

    dep, lamb = bisector(table.lam, norm)
    dep_d, lamb_d = bisector(table.lam, degraded)
    lam0 = table.lam[np.argmin(norm)]
    print(f"最深の線: {lam0:.3f} A, 残留強度 {norm.min():.3f}")
    ok = np.isfinite(lamb)
    if ok.sum() > 2:
        span = (np.nanmax(lamb[ok]) - np.nanmin(lamb[ok])) / lam0 * 2.998e5
        print(f"バイセクタの振れ幅: {span*1000:.0f} m/s")

    if args.save_npz:
        np.savez(args.save_npz, lam=table.lam, mean=mean_spec, norm=norm,
                 degraded=degraded, n=n)
        print(f"wrote {args.save_npz}")

    # --- 図 ------------------------------------------------------------------
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))

    ax[0].plot(table.lam, norm, lw=0.7, label=f"synth R={table.resolving_power:.0e}")
    ax[0].plot(table.lam, degraded, lw=1.4,
               label=f"degraded to R={args.R_inst:.0e}")
    ax[0].set_xlabel(r"$\lambda$ [$\AA$]")
    ax[0].set_ylabel("normalized intensity")
    ax[0].legend(fontsize=8)
    ax[0].set_title("disk-centre intensity (horizontal mean)")

    m = np.abs(table.lam - lam0) < 0.5
    ax[1].plot(table.lam[m], norm[m], lw=1.0)
    ax[1].plot(table.lam[m], degraded[m], lw=1.4)
    ax[1].set_xlabel(r"$\lambda$ [$\AA$]")
    ax[1].set_title(f"line at {lam0:.2f} A")

    c = 2.998e5
    ax[2].plot((lamb - lam0) / lam0 * c * 1000, dep, "-o", ms=3,
               label="synth")
    ax[2].plot((lamb_d - lam0) / lam0 * c * 1000, dep_d, "-s", ms=3,
               label=f"R={args.R_inst:.0e}")
    ax[2].set_xlabel("bisector velocity [m/s]")
    ax[2].set_ylabel("relative depth (0=core, 1=continuum)")
    ax[2].legend(fontsize=8)
    ax[2].set_title("bisector")

    fig.tight_layout()
    fig.savefig(args.out, dpi=140)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
