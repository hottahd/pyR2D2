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


from pyR2D2.synth.core import (OpacityTable, bisector, instrument_profile,
                          synth_column)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--table", required=True, help="高分解能不透明度テーブル (h5)")
    ap.add_argument("--run", default="/scr/a000/c0234hotta/odf/run/d001")
    ap.add_argument("--n", type=int, default=-1, help="出力番号 (-1 で最新)")
    ap.add_argument("--nsnap", type=int, default=1,
                    help="平均するスナップショット数。**p モード (音波) は箱全体で "
                         "同位相なので水平平均では消えず、時間平均でしか落ちない**。"
                         "相関時間 1.5 分に対し出力は 60 秒なので、--snap-step で "
                         "間隔を空けて 10-20 枚averaging するのが目安")
    ap.add_argument("--snap-step", type=int, default=3,
                    help="平均するスナップショットの間隔 (出力番号)")
    ap.add_argument("--stride", type=int, default=4,
                    help="柱を何点おきに使うか (4 なら 64x64 本)")
    ap.add_argument("--xmin", type=float, default=-1.0,
                    help="使う深さの下端 [Mm] (これより深いところは切る)")
    ap.add_argument("--R-inst", type=float, default=3.0e5,
                    help="装置分解能 (岡山の分光器を想定)")
    ap.add_argument("--out", default="notebooks/r2d2_spectrum.png")
    ap.add_argument("--remove-mean-vx", action="store_true",
                    help="各高さで水平平均した vx を引く (箱全体のピストン運動 "
                         "= k_h=0 の成分を落とす)。**k-omega フィルタが使えない "
                         "粗いケーデンスのデータのための代用**。ここで使うのは "
                         "その代用が k-omega とどれだけ一致するかを測るため")
    ap.add_argument("--save-npz", default=None)
    ap.add_argument("--filtered", default=None,
                    help="kw_filter.py の出力 (h5)。指定すると 3D 出力の代わりに "
                         "**音波を落とした**データから合成する。p モードによる "
                         "見かけのシフト (rms 200 m/s) が消えるので、スナップショット "
                         "1 枚でも安定した線位置が出る")
    args = ap.parse_args()

    import pyR2D2

    table = OpacityTable.load(args.table)
    print(f"不透明度テーブル: {table.lam[0]:.2f}-{table.lam[-1]:.2f} A, "
          f"{len(table.lam)} 点, R={table.resolving_power:.3g}, "
          f"vturb={table.meta.get('vturb_kms')} km/s")

    d = pyR2D2.Data(str(Path(args.run) / "data"))
    p = d.p

    # --- フィルタ済みデータを使う場合 ---------------------------------------
    # kw_filter.py は「背景からのずれ」を、高さを間引いた形で保存している。
    # 高さ格子はここに入っている x をそのまま使い、te0 等は p から足す。
    fil = None
    if args.filtered:
        import h5py
        fil = h5py.File(args.filtered, "r")
        x_fil = fil["x"][:]                       # [cm] 昇順 (p.x の部分集合)
        n_fil = fil["n"][:]                       # 元の出力番号
        # p.x のどこに対応するかを引いて te0 等を切り出す
        isel = np.searchsorted(np.array(p.x), x_fil)
        print(f"フィルタ済みデータ: {args.filtered}")
        print(f"  {len(n_fil)} 枚 (n = {n_fil[0]}..{n_fil[-1]}), "
              f"高さ {len(x_fil)} 点, {fil.attrs.get('note', '')}")

    n = int(p.nd) if args.n < 0 else args.n
    print(f"スナップショット: n={n} (t={n*p.dtout/60:.0f} min)")

    snaps = [n - k * args.snap_step for k in range(args.nsnap)]
    snaps = [m for m in snaps if m >= 0][::-1]
    if fil is not None:
        # フィルタ済みファイルに入っている出力番号だけに限る
        snaps = [m for m in snaps if m in set(n_fil.tolist())]
        if not snaps:
            raise SystemExit(f"指定した出力番号がフィルタ済みデータに無い "
                             f"(入っているのは {n_fil[0]}..{n_fil[-1]})")
    print(f"平均するスナップショット: {len(snaps)} 枚 "
          f"(n = {snaps[0]}..{snaps[-1]}, 間隔 {args.snap_step*p.dtout:.0f} s, "
          f"合計 {(snaps[-1]-snaps[0])*p.dtout/60:.0f} 分)")
    x = np.array(p.x) if fil is None else x_fil

    # 上端が先頭になるよう反転し、深いところは切る (tau >> 1 で寄与しない)
    order = np.argsort(x)[::-1]
    x = x[order]
    keep = (x - p.rstar) / 1e8 > args.xmin
    x = x[keep]
    print(f"使う深さ: {len(x)} 点 "
          f"({(x[0]-p.rstar)/1e8:+.2f} .. {(x[-1]-p.rstar)/1e8:+.2f} Mm)")

    t0 = time.time()
    acc = np.zeros(len(table.lam))
    ncol_total = 0
    # **1 枚ごと** のスペクトルも貯める。スナップショット間のばらつきが
    # そのまま「粒状斑起源の RV jitter」の推定になる (docs/12 7 章)。
    per_snap = []
    for si, ns in enumerate(snaps):
        if fil is None:
            d.qr.read(ns, keys=["te", "pr", "ro", "vx"])
            te = np.array(d.qr.te) + np.array(p.te0)[:, None, None]
            pr = np.array(d.qr.pr) + np.array(p.pr0)[:, None, None]
            ro = np.array(d.qr.ro) + np.array(p.ro0)[:, None, None]
            vx = np.array(d.qr.vx)
        else:
            i = int(np.where(n_fil == ns)[0][0])
            te = fil["te"][i] + np.array(p.te0)[isel][:, None, None]
            pr = fil["pr"][i] + np.array(p.pr0)[isel][:, None, None]
            ro = fil["ro"][i] + np.array(p.ro0)[isel][:, None, None]
            vx = fil["vx"][i]
        te, pr, ro, vx = te[order], pr[order], ro[order], vx[order]
        te, pr, ro, vx = te[keep], pr[keep], ro[keep], vx[keep]
        if args.remove_mean_vx:
            vx = vx - vx.mean(axis=(1, 2))[:, None, None]

        ny, nz = te.shape[1], te.shape[2]
        acc_s = np.zeros(len(table.lam))
        ncol_s = 0
        for iy in range(0, ny, args.stride):
            for iz in range(0, nz, args.stride):
                acc_s += synth_column(te[:, iy, iz], pr[:, iy, iz], ro[:, iy, iz],
                                      vx[:, iy, iz], x, table)
                ncol_s += 1
        acc += acc_s
        ncol_total += ncol_s
        per_snap.append(acc_s / ncol_s)
        if si == 0:
            dt = time.time() - t0
            print(f"  1 枚に {dt:.1f} s → 全体 {dt*len(snaps)/60:.1f} 分の見込み")
    mean_spec = acc / ncol_total
    print(f"合成 {time.time()-t0:.1f} s ({ncol_total} 柱 x 波長)")

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
                 degraded=degraded, n=n,
                 per_snap=np.array(per_snap, dtype=np.float32),
                 snaps=np.array(snaps),
                 ncol_per_snap=ncol_total // max(len(snaps), 1),
                 box_mm=float((p.y[-1] - p.y[0]) / 1e8))
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
