"""eval_tau・vertical_upward_rte の試験 (2026-11-06 に直した添字・ゼロ割り・x の倍精度・GIL)。"""

import threading

import numpy as np
import pytest

from pyR2D2.cpp_util import EOS, eval_tau, vertical_upward_rte

SB_OVER_PI = 5.67e-5 / np.pi
TE0 = 5000.0


@pytest.fixture(scope="module")
def eos():
    """不透明度 1、log(te) = log(TE0) + se の表。"""
    log_ro = np.linspace(-30.0, 2.0, 321, dtype=np.float32)
    se = np.linspace(-2.0, 2.0, 401, dtype=np.float32)
    zero = np.zeros((log_ro.size, se.size), dtype=np.float32)
    log_te = (np.log(TE0) + se[None, :] * np.ones_like(zero)).astype(np.float32)
    return EOS(log_ro, se, zero, zero, log_te, zero)


def column(values):
    return np.asarray(values, dtype=np.float32)[:, None, None]


def source(se):
    return SB_OVER_PI * (TE0 * np.exp(se)) ** 4


def test_tau_uniform_has_no_zero_division(eos):
    x = np.arange(16.0)
    ro = column(np.full(16, 1.0e-3))
    tau = eval_tau(ro, np.zeros_like(ro), x, eos)[:, 0, 0]
    assert tau == pytest.approx(1.0e-3 * (x[-1] - x), rel=1.0e-6, abs=1.0e-12)


def test_tau_nearly_equal_neighbours_stay_finite(eos):
    """float32 で 1 ulp だけ違う隣どうし (float の log では差が 0 になり、以前は inf)。"""
    base = np.float32(1.0e-3)
    values = [base, np.nextafter(base, np.float32(1))] * 8
    ro = column(values)
    tau = eval_tau(ro, np.zeros_like(ro), np.arange(16.0), eos)[:, 0, 0]
    assert np.isfinite(tau).all()
    assert tau == pytest.approx(1.0e-3 * (15.0 - np.arange(16.0)), rel=1.0e-5, abs=1.0e-9)


def test_tau_exponential_atmosphere_is_exact(eos):
    x = np.linspace(0.0, 10.0, 41)
    ro = np.exp(-x)
    tau = eval_tau(column(ro), column(np.zeros(41)), x, eos)[:, 0, 0]
    assert tau == pytest.approx(ro - ro[-1], rel=1.0e-5, abs=1.0e-9)


def test_rte_uniform_source_has_no_zero_division(eos):
    ro = column(np.full(16, 1.0e-3))
    rt = vertical_upward_rte(ro, np.zeros_like(ro), np.arange(16.0), eos)
    assert rt[0, 0] == pytest.approx(source(0.0), rel=1.0e-5)


def test_rte_reaches_the_top_face(eos):
    """最上面まで積む: 一番上のセルを変えれば結果が変わり、面の値は解析解に一致する。"""
    n = 40
    x = np.arange(float(n))
    ro = column(np.full(n, 1.0e-1))
    se = np.zeros(n)
    base = vertical_upward_rte(ro, column(se), x, eos)[0, 0]
    se[-1] = 0.5
    changed = vertical_upward_rte(ro, column(se), x, eos)[0, 0]
    # 最上面の S は log の平均で e^{4 * 0.25}。そこへ向けて log S が線形な最後の 1 セルを積む
    dtu = 1.0e-1
    s0, s1 = source(0.0), source(0.25)
    expected = s0 * np.exp(-dtu) + dtu * (s1 - s0 * np.exp(-dtu)) / (np.log(s1 / s0) + dtu)
    assert base == pytest.approx(s0, rel=1.0e-5)
    assert changed == pytest.approx(expected, rel=1.0e-5)


def test_rte_source_term_denominator_zero(eos):
    """log S の変化が dtu をちょうど打ち消すとき (分母 0)。解は S0 e^{-t} (1 + t)。"""
    n = 41
    ro_value, dx = 1.0e-2, 5.0
    dtu = ro_value * dx
    se = 1.0 - 0.25 * dtu * np.arange(n)  # log S = 4 se なので 1 セルで -dtu
    rt = vertical_upward_rte(column(np.full(n, ro_value)), column(se), dx * np.arange(float(n)), eos)[0, 0]
    t = dtu * (n - 2)
    s0 = source(0.5 * (se[0] + se[1]))
    assert np.isfinite(rt)
    assert rt == pytest.approx(s0 * np.exp(-t) * (1.0 + t), rel=1.0e-3)


def test_thick_cells_stay_finite(eos):
    """1 セルの厚みが大きくても (e^{-dtu} = 0)、NaN にならず I = S に近づく。"""
    n = 8
    ro = column(np.full(n, 1.0))
    se = np.linspace(0.2, 0.0, n)
    rt = vertical_upward_rte(ro, column(se), 1.0e4 * np.arange(float(n)), eos)[0, 0]
    assert rt == pytest.approx(source(0.5 * (se[-2] + se[-1])), rel=1.0e-2)
    assert np.isfinite(eval_tau(ro, column(se), 1.0e4 * np.arange(float(n)), eos)).all()


def test_x_can_be_float64_with_a_large_offset(eos):
    """float32 なら刻みが格子幅の 0.7% になる大きな座標でも、差分は倍精度で取る。"""
    n = 64
    dx = 2.4e6
    x = dx * np.arange(float(n))
    rng = np.random.default_rng(0)
    ro = column(np.exp(-x / 3.0e7) * (1.0 + 0.01 * rng.random(n)) * 1.0e-7)
    se = column(0.1 * rng.random(n))
    offset = 1.35e11
    # float32 の座標なら 1e-3 程度ずれる。丸めの範囲 (1 ulp 程度) で一致すること
    assert eval_tau(ro, se, x + offset, eos) == pytest.approx(eval_tau(ro, se, x, eos), rel=1.0e-6)
    assert vertical_upward_rte(ro, se, x + offset, eos) == pytest.approx(vertical_upward_rte(ro, se, x, eos), rel=1.0e-6)
    # float32 も今までどおり渡せる
    assert eval_tau(ro, se, x.astype(np.float32), eos).shape == ro.shape


def test_too_short_grid_is_an_error(eos):
    ro = column([1.0e-3])
    with pytest.raises(RuntimeError):
        vertical_upward_rte(ro, np.zeros_like(ro), np.zeros(1), eos)
    with pytest.raises(RuntimeError):
        eval_tau(column([1.0e-3, 1.0e-3]), column([0.0, 0.0]), np.zeros(3), eos)


@pytest.mark.parametrize(
    "call",
    [
        lambda ro, se, x, eos: eval_tau(ro, se, x, eos),
        lambda ro, se, x, eos: vertical_upward_rte(ro, se, x, eos),
        lambda ro, se, x, eos: eos.eval(ro, se, "op"),
    ],
    ids=["eval_tau", "vertical_upward_rte", "eos_eval"],
)
def test_same_result_from_threads(eos, call):
    """GIL を放す: 同じ EOS を別スレッドから同時に使っても結果が同じ。"""
    rng = np.random.default_rng(1)
    shape = (64, 32, 32)
    ro = (1.0e-7 * (1.0 + rng.random(shape))).astype(np.float32)
    se = (0.2 * rng.random(shape)).astype(np.float32)
    x = np.arange(64.0)
    expected = call(ro, se, x, eos)
    out = [None] * 8

    def work(i):
        out[i] = call(ro, se, x, eos)

    threads = [threading.Thread(target=work, args=(i,)) for i in range(8)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    assert all(np.array_equal(o, expected) for o in out)
