"""pyR2D2.synth の試験（2026-10-02 に ODF-radiation の tests/test_physics.py から、合成の 11 件と補助をそのまま移した。R2D2plus DEC-589）。"""
import numpy as np
import pytest

pytest.importorskip("h5py")


def test_formal_solution_uses_logarithmic_mean():
    """log スキームの光学的厚みが**対数平均**になっていること。

    alpha が 1 セルで 1 -> 100 と変わる場合:
        算術平均 50.5 / 幾何平均 10 / **対数平均 21.5**
    ln(alpha) が経路長について線形 (指数関数的変化) なら 21.5 が厳密値で、
    幾何平均の 10 は「log 空間の中点を代表値にする」別物 (過小評価)。
    """
    import numpy as np

    from pyR2D2.synth import formal_solution

    a0, a1, ds = 1.0, 100.0, 1.0
    exact = (a1 - a0) / np.log(a1 / a0)          # 対数平均 = 21.497
    s = np.linspace(0.0, ds, 200001)
    num = np.trapezoid(a0 * np.exp(s / ds * np.log(a1 / a0)), s) / ds
    assert abs(num - exact) / exact < 1e-6

    # 2 点 1 セルの解析解と突き合わせる。深い側 (添字 1) が上流で、
    # そこは拡散極限として I = S_u で初期化される。
    #     I = S_u exp(-d) + (S_d - S_u exp(-d)) / (1 + ln(S_d/S_u)/d)
    ds2, s_u, s_d = 0.01, 2.0, 1.0
    alpha = np.array([[a0], [a1]])
    source = np.array([[s_d], [s_u]])          # 上端が先頭
    I = formal_solution(alpha, source, np.array([ds2]), scheme="log")[0]

    def analytic(d):
        ex = np.exp(-d)
        return s_u * ex + (s_d - s_u * ex) / (1.0 + np.log(s_d / s_u) / d)

    d_log = exact * ds2                        # 対数平均を使った場合
    d_lin = 0.5 * (a0 + a1) * ds2              # 算術平均を使った場合
    assert abs(I - analytic(d_log)) / I < 1e-10
    # 算術平均だと有意に違う値になる (= 取り違えたら気付ける)
    assert abs(I - analytic(d_lin)) / I > 1e-2


def test_formal_solution_schemes_agree_on_fine_grid():
    """細かい格子では linear と log が同じ答えに収束すること。"""
    import numpy as np

    from pyR2D2.synth import formal_solution

    n = 4000
    x = np.linspace(0.0, 1.0, n)
    # 光球を模した指数関数的な alpha と、緩やかに変わる S
    alpha = np.exp(8.0 * x)[:, None]
    source = (1.0 + 3.0 * x)[:, None]
    dx = np.diff(x)
    I_log = formal_solution(alpha[::-1], source[::-1], dx, scheme="log")[0]
    I_lin = formal_solution(alpha[::-1], source[::-1], dx, scheme="linear")[0]
    assert abs(I_log - I_lin) / I_log < 1.0e-3


def test_formal_solution_rejects_negative_input():
    """負の alpha / S を log スキームに渡したら黙って NaN を返さず落ちること。"""
    import numpy as np
    import pytest

    from pyR2D2.synth import formal_solution

    alpha = np.array([[1.0], [-1.0]])
    source = np.array([[1.0], [1.0]])
    with pytest.raises(ValueError, match="正の alpha"):
        formal_solution(alpha, source, np.array([1.0]), scheme="log")


def test_doppler_shift_direction():
    """観測者へ向かう速度が**青方偏移**になること。

    2026-08-17 まで符号が逆で、上昇流が赤方偏移していた。その結果
    対流ブルーシフトが符号ごと反転し、合成スペクトルが観測より
    数百 m/s 赤いという誤った結論を出していた。
    """
    import numpy as np

    from pyR2D2.synth import doppler_shift_rows

    R = 1.0e6
    nlam = 2001
    # 中央に 1 本だけ深い線がある不透明度
    logk = np.full((1, nlam), -3.0)
    logk[0, nlam // 2] = 1.0
    c = 2.99792458e10

    for v_kms, sign in ((+3.0, -1), (-3.0, +1)):
        out = doppler_shift_rows(logk, np.array([v_kms * 1e5]), R)
        i = int(np.argmax(out[0]))
        moved = i - nlam // 2                  # 正 = 長波長 (赤) へ移動
        assert np.sign(moved) == sign, (
            f"v = {v_kms:+.1f} km/s で線が "
            f"{'赤' if moved > 0 else '青'}へ動いた (期待は逆)")
        expect = abs(np.log1p(v_kms * 1e5 / c) / np.log1p(1.0 / R))
        assert abs(abs(moved) - expect) <= 1.0


def _toy_table(nlam=64):
    """試験用の小さな不透明度テーブル。線を 1 本入れておく。"""
    import numpy as np

    from pyR2D2.synth import OpacityTable

    T = np.logspace(np.log10(3500.0), np.log10(9000.0), 8)
    P = np.logspace(2.0, 6.0, 6)
    lam = 6300.0 * (1.0 + 1.0 / 3.0e5) ** np.arange(nlam)
    logk = np.zeros((len(T), len(P), nlam))
    for i, t in enumerate(T):
        for j, p in enumerate(P):
            # 連続吸収は T, P に滑らかに依存、中央に線を 1 本
            logk[i, j] = -2.0 + 0.5 * np.log10(p / 1e4) + np.log10(t / 5000.0)
            logk[i, j, nlam // 2] += 2.0
            logk[i, j, nlam // 2 - 1] += 1.0
            logk[i, j, nlam // 2 + 1] += 1.0
    return OpacityTable(T=T, P=P, lam=lam, logk=logk, meta={})


def _toy_atmosphere(nx=24, ny=8, nz=8, horizontal=True, seed=3):
    """試験用の 3D 大気。`horizontal=False` なら水平一様。"""
    import numpy as np

    rng = np.random.default_rng(seed)
    # 上端が先頭 (降順)
    x = np.linspace(0.5e8, -0.5e8, nx)
    z0 = np.linspace(0.0, 1.0, nx)
    T1 = 4200.0 + 4200.0 * z0
    P1 = 10.0 ** (2.5 + 3.0 * z0)
    R1 = 10.0 ** (-9.5 + 2.5 * z0)
    T3 = np.repeat(T1[:, None, None], ny, 1).repeat(nz, 2).copy()
    P3 = np.repeat(P1[:, None, None], ny, 1).repeat(nz, 2).copy()
    R3 = np.repeat(R1[:, None, None], ny, 1).repeat(nz, 2).copy()
    v = [np.zeros_like(T3) for _ in range(3)]
    if horizontal:
        # 水平方向に滑らかな揺らぎ (周期的であること)
        jj = np.arange(ny)[None, :, None]
        kk = np.arange(nz)[None, None, :]
        f = (np.sin(2 * np.pi * jj / ny) * np.cos(2 * np.pi * kk / nz))
        T3 *= 1.0 + 0.08 * f
        P3 *= 1.0 + 0.05 * f
        R3 *= 1.0 + 0.05 * f
        v[0] = 1.0e5 * f * np.ones_like(T3)
        v[1] = 0.5e5 * np.cos(2 * np.pi * jj / ny) * np.ones_like(T3)
        v[2] = 0.3e5 * np.sin(2 * np.pi * kk / nz) * np.ones_like(T3)
    return T3, P3, R3, v[0], v[1], v[2], x


def test_synth_ray_reduces_to_synth_column_at_mu_one():
    """**mu=1 では既存の synth_column と一致すること** (最重要の回帰試験)。

    傾いた光線の実装が、垂直入射の極限で従来の柱の計算に戻らなければ
    どこかが間違っている。
    """
    import numpy as np

    from pyR2D2.synth import synth_column, synth_ray

    tab = _toy_table()
    T3, P3, R3, vx, vy, vz, x = _toy_atmosphere()
    dy = dz = 4.0e7
    for (jy, jz) in ((0, 0), (3, 5)):
        I_ray = synth_ray(T3, P3, R3, vx, vy, vz, x, dy, dz, tab,
                          mu=1.0, phi=0.0, iy0=jy, iz0=jz)
        I_col = synth_column(T3[:, jy, jz], P3[:, jy, jz], R3[:, jy, jz],
                             vx[:, jy, jz], x, tab)
        assert np.allclose(I_ray, I_col, rtol=1e-12), (
            f"mu=1 で一致しない (柱 {jy},{jz}): "
            f"最大相対差 {np.max(np.abs(I_ray/I_col - 1)):.2e}")


def test_synth_ray_matches_plane_parallel_when_horizontally_uniform():
    """**水平一様な大気なら、傾けても平行平面の答えと一致すること**。

    水平構造が無ければ「別の柱を通る」効果は消え、経路長が 1/mu 倍に
    なるだけになる。これは `formal_solution(mu=...)` の平行平面の扱いと
    一致しなければならない。**幾何 (経路長) が正しいことの検証**である。
    """
    import numpy as np

    from pyR2D2.synth import synth_column, synth_ray

    tab = _toy_table()
    T3, P3, R3, vx, vy, vz, x = _toy_atmosphere(horizontal=False)
    dy = dz = 4.0e7
    for mu in (0.9, 0.6, 0.3):
        I_ray = synth_ray(T3, P3, R3, vx, vy, vz, x, dy, dz, tab,
                          mu=mu, phi=0.7, iy0=2, iz0=1)
        I_pp = synth_column(T3[:, 2, 1], P3[:, 2, 1], R3[:, 2, 1],
                            vx[:, 2, 1], x, tab, mu=mu)
        assert np.allclose(I_ray, I_pp, rtol=1e-10), (
            f"mu={mu} で平行平面と一致しない: "
            f"最大相対差 {np.max(np.abs(I_ray/I_pp - 1)):.2e}")


def test_synth_ray_is_periodic():
    """**開始位置を箱 1 周ぶんずらしても同じ答えになること**。

    水平の巻き (剰余) が正しく効いていることの検証。`fmod` は負の値で
    負を返すので、素朴に書くと添字が範囲外になる。
    """
    import numpy as np

    from pyR2D2.synth import synth_ray

    tab = _toy_table()
    T3, P3, R3, vx, vy, vz, x = _toy_atmosphere()
    ny, nz = T3.shape[1], T3.shape[2]
    dy = dz = 4.0e7
    kw = dict(mu=0.4, phi=1.1, scheme="log")
    I0 = synth_ray(T3, P3, R3, vx, vy, vz, x, dy, dz, tab, iy0=1.5, iz0=2.5, **kw)
    for shift in (+1, -2):
        I1 = synth_ray(T3, P3, R3, vx, vy, vz, x, dy, dz, tab,
                       iy0=1.5 + shift * ny, iz0=2.5 - shift * nz, **kw)
        assert np.allclose(I0, I1, rtol=1e-12), (
            f"箱 {shift} 周ずらすと答えが変わる: "
            f"最大相対差 {np.max(np.abs(I1/I0 - 1)):.2e}")


def test_substeps_needed():
    """細分数が「水平に 1 格子ずれる前に刻む」条件を満たすこと。"""
    import numpy as np

    from pyR2D2.synth import substeps_needed

    x = np.linspace(0.5e8, -0.5e8, 26)      # 鉛直 40 km 格子
    dx = abs(x[1] - x[0])
    dy = dx                                  # 水平も 40 km (現実的な比)
    assert substeps_needed(x, 1.0, dy, dy) == 1
    seen_gt1 = False
    for mu in (0.9, 0.7, 0.5, 0.3, 0.1):
        n = substeps_needed(x, mu, dy, dy)
        st = np.sqrt(1 - mu ** 2)
        assert dx * st / mu / n <= dy * 1.0000001, (
            f"mu={mu}: {n} 分割しても水平に {dx*st/mu/n/1e5:.0f} km ずれる "
            f"(格子 {dy/1e5:.0f} km)")
        if n > 1:
            seen_gt1 = True
    # **試験が実質的に効いていることの確認**: どこかで 2 分割以上が要るはず。
    # (最初 dy を dx の 10 倍に書いてしまい、全部 n=1 で素通りしていた)
    assert seen_gt1, "全ての mu で n=1 になった。試験の格子設定が甘い"


def test_synth_ray_velocity_weights():
    """**視線速度の重み `wvert`/`whoriz` が意図どおり効くこと**。

    限界効果の分解 (`docs/12` 13.8 節) では「水平速度を消した合成」を
    作る。既定値では従来と 1 ビットも変わらないこと、0 にすると本当に
    その成分が消えること、そして **mu<1 では実際に答えが変わること**
    (試験が素通りしていないこと) を確かめる。
    """
    import numpy as np

    from pyR2D2.synth import synth_ray

    tab = _toy_table()
    T3, P3, R3, vx, vy, vz, x = _toy_atmosphere()
    dy = dz = 4.0e7
    zero = np.zeros_like(vx)
    kw = dict(phi=0.7, iy0=2, iz0=1)

    # 1. 既定値は従来と同一
    I_def = synth_ray(T3, P3, R3, vx, vy, vz, x, dy, dz, tab, mu=0.5, **kw)
    I_one = synth_ray(T3, P3, R3, vx, vy, vz, x, dy, dz, tab, mu=0.5,
                      wvert=1.0, whoriz=1.0, **kw)
    assert np.array_equal(I_def, I_one), "既定値で答えが変わった"

    # 2. 両方 0 なら速度をゼロにしたのと同じ
    I_off = synth_ray(T3, P3, R3, vx, vy, vz, x, dy, dz, tab, mu=0.5,
                      wvert=0.0, whoriz=0.0, **kw)
    I_novel = synth_ray(T3, P3, R3, zero, zero, zero, x, dy, dz, tab,
                        mu=0.5, **kw)
    assert np.allclose(I_off, I_novel, rtol=1e-12), "重み 0 で速度が消えていない"

    # 3. mu=1 では sqrt(1-mu^2)=0 なので水平成分は元々効かない
    a = synth_ray(T3, P3, R3, vx, vy, vz, x, dy, dz, tab, mu=1.0, **kw)
    b = synth_ray(T3, P3, R3, vx, vy, vz, x, dy, dz, tab, mu=1.0,
                  whoriz=0.0, **kw)
    assert np.allclose(a, b, rtol=1e-12), "mu=1 で水平成分が効いてしまった"

    # 4. **mu<1 では水平成分を消すと答えが変わること** (試験が効いている確認)
    c = synth_ray(T3, P3, R3, vx, vy, vz, x, dy, dz, tab, mu=0.4,
                  whoriz=0.0, **kw)
    d = synth_ray(T3, P3, R3, vx, vy, vz, x, dy, dz, tab, mu=0.4, **kw)
    assert np.max(np.abs(c / d - 1)) > 1e-4, (
        "mu=0.4 で水平成分を消しても答えが変わらない。"
        "重みが繋がっていない疑い")


def test_ray_contribution_splits_v_los():
    """`ray_contribution` が視線速度を鉛直・水平に正しく分けること。"""
    import numpy as np

    from pyR2D2.synth import ray_contribution

    tab = _toy_table()
    T3, P3, R3, _, _, _, x = _toy_atmosphere()
    dy = dz = 4.0e7
    A, B, C = 1.3e5, -0.7e5, 0.4e5          # 一様な速度場 [cm/s]
    vx = np.full_like(T3, A)
    vy = np.full_like(T3, B)
    vz = np.full_like(T3, C)
    mu, phi = 0.4, 0.9
    st = np.sqrt(1.0 - mu ** 2)
    r = ray_contribution(T3, P3, R3, vx, vy, vz, x, dy, dz, tab,
                         [len(tab.lam) // 2], mu=mu, phi=phi, iy0=1, iz0=2)
    assert np.allclose(r["v_vert"], mu * A)
    assert np.allclose(r["v_horiz"], st * (B * np.cos(phi) + C * np.sin(phi)))
    # tau は上端 0 から単調増加
    assert r["tau"][0, 0] == 0.0
    assert np.all(np.diff(r["tau"][:, 0]) > 0.0)


def test_ray_contribution_integrates_to_intensity():
    """**寄与関数を経路長で積むと出射強度になること**。

        I = int CF ds + S(下端) exp(-tau(下端))

    第 2 項は形式解が拡散極限で置く下端の初期値。これを足さないと
    「光学的に薄い試験大気で 3 割足りない」ことになり、重みとしての
    妥当性が判断できない。細かい格子で数 % に収まればよい
    (log スキームの形式解と台形積分は厳密には一致しない)。
    """
    import numpy as np

    from pyR2D2.synth import ray_contribution, synth_ray

    tab = _toy_table()
    T3, P3, R3, vx, vy, vz, x = _toy_atmosphere(nx=200)
    dy = dz = 4.0e7
    idx = [len(tab.lam) // 2, 2]            # 線コアと連続光
    for mu in (1.0, 0.6):
        r = ray_contribution(T3, P3, R3, vx, vy, vz, x, dy, dz, tab, idx,
                             mu=mu, phi=0.3, iy0=1, iz0=2)
        I = synth_ray(T3, P3, R3, vx, vy, vz, x, dy, dz, tab,
                      mu=mu, phi=0.3, iy0=1, iz0=2)[idx]
        # 下端の初期値 (拡散極限) を足す
        from pyR2D2.synth._physics import planck_lambda
        # **傾いた光線の下端は別の柱にあるので、r["T"] の下端を使うこと**
        # (T3[-1,1,2] を使うと mu=0.6 で 15% 外れる)
        S_bot = planck_lambda(tab.lam[idx], float(r["T"][-1]))
        I_cf = (r["cf"] * r["w"][:, None]).sum(axis=0) \
            + S_bot * np.exp(-r["tau"][-1])
        assert np.allclose(I_cf, I, rtol=0.03), (
            f"mu={mu}: CF の積分が強度と合わない "
            f"{I_cf} vs {I} (比 {I_cf/I})")
