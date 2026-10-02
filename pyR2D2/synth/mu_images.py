#!/usr/bin/env python3
"""**傾けた視線**での出射強度の画像を作り、周縁減光を観測と比べる。

**なぜこれが検証になるか**: 傾いた光線の実装 (`pyR2D2.synth.synth_ray`) が
正しければ、mu を振って得た強度から**周縁減光**が出る。そして

    <I> / I(mu=1) = 2 * int_0^1 I(mu) mu dmu / I(1)

は Neckel の Hamburg アトラスの**絶対値**から直接読める
(全面平均 file01-07 と視線中心 file11-17 の連続光の比)。
6300 A で **0.828** (実測)。**これは合成の他の部分を一切使わない独立な検定**
である。

出力:
    * mu ごとの連続光強度の画像 (粒状斑が斜めから見えるとどうなるか)
    * 強度のコントラスト (rms/平均) の mu 依存
    * 周縁減光の曲線と、アトラスから読んだ値との比較

使い方:
    python scripts/synth_mu_images.py --table $W/opac_fe6302_vald.h5 \
        --n 141 --stride 2 --out notebooks/mu_images.png
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np


from pyR2D2.synth.core import OpacityTable, substeps_needed, synth_ray

MUS = (1.0, 0.8, 0.6, 0.4, 0.2)


def continuum_table(table, lam_lo, lam_hi):
    """連続光の狭い窓だけに切り詰めたテーブルを返す (計算量を落とすため)。"""
    m = (table.lam >= lam_lo) & (table.lam <= lam_hi)
    return OpacityTable(T=table.T, P=table.P, lam=table.lam[m],
                        logk=table.logk[:, :, m], meta=table.meta)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--table", required=True)
    ap.add_argument("--run", default="/scr/a000/c0234hotta/odf/run/d001")
    ap.add_argument("--n", type=int, default=141)
    ap.add_argument("--stride", type=int, default=2,
                    help="光線を何格子おきに引くか")
    ap.add_argument("--lam-lo", type=float, default=6285.4)
    ap.add_argument("--lam-hi", type=float, default=6285.7)
    ap.add_argument("--xmin", type=float, default=-0.5)
    ap.add_argument("--nphi", type=int, default=4,
                    help="周縁減光を測るときに平均する方位角の数")
    ap.add_argument("--out", default="notebooks/mu_images.png")
    args = ap.parse_args()

    import pyR2D2

    tab = continuum_table(OpacityTable.load(args.table), args.lam_lo, args.lam_hi)
    print(f"連続光の窓: {tab.lam[0]:.3f}-{tab.lam[-1]:.3f} A ({len(tab.lam)} 点)")

    d = pyR2D2.Data(str(Path(args.run) / "data"))
    p = d.p
    d.qr.read(args.n, keys=["te", "pr", "ro", "vx", "vy", "vz"])
    te = np.array(d.qr.te) + np.array(p.te0)[:, None, None]
    pr = np.array(d.qr.pr) + np.array(p.pr0)[:, None, None]
    ro = np.array(d.qr.ro) + np.array(p.ro0)[:, None, None]
    vx, vy, vz = (np.array(getattr(d.qr, k)) for k in ("vx", "vy", "vz"))
    x = np.array(p.x)
    o = np.argsort(x)[::-1]
    x = x[o]
    te, pr, ro, vx, vy, vz = (a[o] for a in (te, pr, ro, vx, vy, vz))
    keep = (x - p.rstar) / 1e8 > args.xmin
    x = x[keep]
    te, pr, ro, vx, vy, vz = (a[keep] for a in (te, pr, ro, vx, vy, vz))
    dy = float(p.y[1] - p.y[0])
    dz = float(p.z[1] - p.z[0])
    ny, nz = te.shape[1], te.shape[2]
    print(f"格子 {len(x)} x {ny} x {nz}, 水平 {dy/1e5:.0f} km, "
          f"箱 {ny*dy/1e8:.2f} Mm")

    ys = np.arange(0, ny, args.stride)
    zs = np.arange(0, nz, args.stride)
    imgs, means, contrasts, subs = {}, {}, {}, {}
    for mu in MUS:
        nsub = substeps_needed(x, mu, dy, dz)
        subs[mu] = nsub
        t0 = time.time()
        img = np.empty((len(ys), len(zs)))
        for a, iy in enumerate(ys):
            for b, iz in enumerate(zs):
                I = synth_ray(te, pr, ro, vx, vy, vz, x, dy, dz, tab,
                              mu=mu, phi=0.0, iy0=iy, iz0=iz, nsub=nsub)
                img[a, b] = I.mean()
        imgs[mu] = img
        means[mu] = img.mean()
        contrasts[mu] = img.std() / img.mean()
        print(f"  mu={mu:.1f} (細分 {nsub}): 平均 {means[mu]:.4e}, "
              f"コントラスト {100*contrasts[mu]:.2f}%  ({time.time()-t0:.0f} s)")

    # --- 周縁減光の検定 -----------------------------------------------------
    # 方位角を平均する (1 方向だけだと箱の中の流れの偏りが残る)
    mu_fine = np.array([1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1])
    Imu = np.empty(len(mu_fine))
    print("\n周縁減光を測る (方位角 %d 点で平均):" % args.nphi)
    for i, mu in enumerate(mu_fine):
        nsub = substeps_needed(x, mu, dy, dz)
        acc, nc = 0.0, 0
        for ph in np.linspace(0.0, 2 * np.pi, args.nphi, endpoint=False):
            for iy in range(0, ny, 4 * args.stride):
                for iz in range(0, nz, 4 * args.stride):
                    acc += synth_ray(te, pr, ro, vx, vy, vz, x, dy, dz, tab,
                                     mu=mu, phi=ph, iy0=iy, iz0=iz,
                                     nsub=nsub).mean()
                    nc += 1
        Imu[i] = acc / nc
        print(f"  mu={mu:.1f}: I/I(1) = {Imu[i]/Imu[0]:.4f}")

    # <I>/I(1) = 2 int I(mu) mu dmu / I(1)。mu=0 は I=0 とみなす
    mm = np.concatenate([[0.0], mu_fine[::-1]])
    II = np.concatenate([[0.0], Imu[::-1]])
    ratio = 2.0 * np.trapezoid(II * mm, mm) / Imu[0]
    print(f"\n**<I>/I(mu=1) = {ratio:.3f}**")
    print("  アトラスの絶対値から読んだ観測値 = **0.828** "
          "(6300 A、全面平均 0.2535 / 視線中心 0.3061 W/cm2/sr/A)")
    print(f"  差 {100*(ratio/0.828-1):+.1f}%")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(16, 7.6))
    gs = fig.add_gridspec(2, len(MUS), height_ratios=[1.15, 1], hspace=0.42,
                          wspace=0.28)
    ext = [0, nz * dz / 1e8, 0, ny * dy / 1e8]
    for i, mu in enumerate(MUS):
        ax = fig.add_subplot(gs[0, i])
        # **各パネルを自分の平均で規格化する**。全部を I(mu=1) で割ると
        # 低 mu が真っ黒になって構造が見えない (最初それで失敗した)。
        im = imgs[mu] / means[mu]
        ax.imshow(im, origin="lower", extent=ext, cmap="afmhot",
                  vmin=0.62, vmax=1.38)
        ax.set_title(f"mu = {mu:.1f}\nI/I(1) = {means[mu]/means[1.0]:.3f}, "
                     f"rms {100*contrasts[mu]:.1f}%", fontsize=10)
        ax.set_xlabel("z [Mm]", fontsize=9)
        if i == 0:
            ax.set_ylabel("y [Mm]", fontsize=9)
        ax.tick_params(labelsize=8)

    np.savez(str(Path(args.out).with_suffix(".npz")),
             mus=np.array(MUS), mu_fine=mu_fine, Imu=Imu, ratio=ratio,
             contrasts=np.array([contrasts[m] for m in MUS]),
             means=np.array([means[m] for m in MUS]),
             **{f"img{m:.1f}": imgs[m] for m in MUS})

    ax = fig.add_subplot(gs[1, :2])
    ax.plot(mu_fine, Imu / Imu[0], "o-", label="synth")
    ax.set_xlabel("mu"); ax.set_ylabel("I(mu) / I(1)")
    ax.set_title(f"limb darkening   <I>/I(1) = {ratio:.3f} "
                 f"(obs 0.828)")
    ax.grid(alpha=0.3); ax.legend(fontsize=9)

    ax = fig.add_subplot(gs[1, 3:])
    ax.plot(list(MUS), [100 * contrasts[m] for m in MUS], "s-")
    ax.set_xlabel("mu"); ax.set_ylabel("intensity contrast [%]")
    ax.set_title("granulation contrast vs mu")
    ax.grid(alpha=0.3)

    out = Path(args.out); out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=120)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
