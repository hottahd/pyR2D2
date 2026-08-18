#!/usr/bin/env python3
"""**円盤積分** — 太陽を「恒星として」見たスペクトルを合成する。

**なぜこれが本命か**: 佐藤さんの分光器が見るのは円盤全面の光である。
mu=1 の柱だけを並べても観測量にはならない。円盤の半分の面積は
mu < 0.71 にあり、そこでは水平速度が視線に乗り、自転も効く。

**幾何**: 円盤上の点を (mu, psi) で表す。mu = cos(視線と局所鉛直の角)、
psi は円盤中心まわりの方位角。投影面積要素は

    dA_proj = R^2 mu dmu dpsi     (積分すると pi R^2 で正しい)

したがって

    F(lambda) ∝ int_0^2pi int_0^1 I(mu, psi; lambda) mu dmu dpsi

**自転**: 自転軸が視線に垂直 (太陽は 7 度傾いているだけ) とすると、
円盤上の点の視線速度は

    v_rot(mu, psi) = v_eq sqrt(1-mu^2) sin(psi)

太陽の赤道自転速度は約 **2.0 km/s** (6300 A で 42 mA)。線の半値全幅が
150 mA 程度なので**無視できない**。

**方位角の扱い**: 箱は水平方向に統計的に等方なので、ある mu に対する
スペクトルを局所方位角 phi で平均すれば、円盤方位角 psi によらない。
よって「phi で平均した I(mu)」を作ってから、psi 積分は**自転シフトの
ためだけ**に回せばよい。これで計算量が (mu x phi) だけで済む。

**検証**: 出来上がりは Neckel アトラスの**全面平均** (file01-07) と
比べる。2026-08-17 まで mu=1 の合成をこれと比べていたのは誤りだったが、
**円盤積分にとってはこちらが正しい相手**である。

使い方:
    python scripts/synth_disk.py --table $W/opac_fe6302_vald.h5 \
        --n 141 --nmu 6 --nphi 4 --stride 16 --save-npz $W/disk.npz
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from odfgen.synth import OpacityTable, substeps_needed, synth_ray

C_KMS = 2.99792458e5


def shift_spectrum(lam, spec, v_kms):
    """等分解能格子の上でスペクトルを視線速度 [km/s] だけずらす。

    観測者へ向かう (v>0) と青方偏移。`doppler_shift_rows` と同じ約束。
    """
    if v_kms == 0.0:
        return spec
    dln = np.log(lam[1] / lam[0])
    sh = np.log1p(v_kms / C_KMS) / dln
    idx = np.arange(len(lam))
    return np.interp(idx + sh, idx, spec, left=spec[0], right=spec[-1])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--table", required=True)
    ap.add_argument("--run", default="/scr/a000/c0234hotta/odf/run/d001")
    ap.add_argument("--n", type=int, default=141)
    ap.add_argument("--nsnap", type=int, default=1)
    ap.add_argument("--snap-step", type=int, default=3)
    ap.add_argument("--nmu", type=int, default=6,
                    help="mu 方向の Gauss 求積の点数。**重み mu dmu 込みで "
                         "積分するので 6 点でも精度が出る**")
    ap.add_argument("--nphi", type=int, default=4,
                    help="局所方位角の点数 (箱の中の流れの偏りを均す)")
    ap.add_argument("--stride", type=int, default=16,
                    help="各 (mu, phi) で引く光線の間隔 [格子]")
    ap.add_argument("--npsi", type=int, default=24,
                    help="自転を入れるときの円盤方位角の点数")
    ap.add_argument("--vrot", type=float, default=2.0,
                    help="赤道自転速度 [km/s]。0 にすると自転を入れない")
    ap.add_argument("--xmin", type=float, default=-0.5)
    ap.add_argument("--save-npz", default=None)
    args = ap.parse_args()

    import pyR2D2

    tab = OpacityTable.load(args.table)
    print(f"テーブル {tab.lam[0]:.1f}-{tab.lam[-1]:.1f} A, {len(tab.lam)} 点")

    d = pyR2D2.Data(str(Path(args.run) / "data"))
    p = d.p
    dy = float(p.y[1] - p.y[0])
    dz = float(p.z[1] - p.z[0])

    # mu の求積: int_0^1 f(mu) mu dmu を Gauss-Legendre で
    gx, gw = np.polynomial.legendre.leggauss(args.nmu)
    mus = 0.5 * (gx + 1.0)                 # [0,1] へ
    wmu = 0.5 * gw * mus                   # 重み mu dmu 込み
    wmu = wmu / wmu.sum()
    print(f"mu の求積点: {np.round(mus,3)}")
    print(f"  重み (mu dmu 込み、規格化): {np.round(wmu,3)}")

    snaps = [args.n - k * args.snap_step for k in range(args.nsnap)][::-1]
    Imu = np.zeros((args.nmu, len(tab.lam)))
    t0 = time.time()
    for ns in snaps:
        d.qr.read(ns, keys=["te", "pr", "ro", "vx", "vy", "vz"])
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
        ny, nz = te.shape[1], te.shape[2]

        for im, mu in enumerate(mus):
            nsub = substeps_needed(x, mu, dy, dz)
            acc = np.zeros(len(tab.lam))
            nc = 0
            for ph in np.linspace(0.0, 2 * np.pi, args.nphi, endpoint=False):
                for iy in range(0, ny, args.stride):
                    for iz in range(0, nz, args.stride):
                        acc += synth_ray(te, pr, ro, vx, vy, vz, x, dy, dz,
                                         tab, mu=mu, phi=ph, iy0=iy, iz0=iz,
                                         nsub=nsub)
                        nc += 1
            Imu[im] += acc / nc
            print(f"  n={ns} mu={mu:.3f} (細分 {nsub}, {nc} 本) "
                  f"{time.time()-t0:.0f} s", flush=True)
    Imu /= len(snaps)

    # --- 円盤積分 -----------------------------------------------------------
    if args.vrot > 0:
        psis = np.linspace(0.0, 2 * np.pi, args.npsi, endpoint=False)
        F = np.zeros(len(tab.lam))
        for im, mu in enumerate(mus):
            st = np.sqrt(max(0.0, 1.0 - mu * mu))
            for ps in psis:
                v = args.vrot * st * np.sin(ps)
                F += wmu[im] / args.npsi * shift_spectrum(tab.lam, Imu[im], v)
        print(f"\n自転 {args.vrot} km/s を入れた "
              f"(円盤方位角 {args.npsi} 点)")
    else:
        F = (wmu[:, None] * Imu).sum(axis=0)
        print("\n自転なし")

    norm = F / np.percentile(F, 99)
    I1 = Imu[np.argmax(mus)]
    print(f"円盤積分/視線中心 の連続光比 = "
          f"{np.percentile(F,99)/np.percentile(I1,99):.3f}  "
          f"(mu={mus.max():.3f} を視線中心の代用とした概算)")
    print(f"最深の線: {tab.lam[np.argmin(norm)]:.3f} A, "
          f"残留強度 {norm.min():.3f}")

    if args.save_npz:
        np.savez(args.save_npz, lam=tab.lam, mean=F, norm=norm,
                 degraded=norm, n=args.n, mus=mus, wmu=wmu, Imu=Imu,
                 vrot=args.vrot)
        print(f"wrote {args.save_npz}")


if __name__ == "__main__":
    main()
