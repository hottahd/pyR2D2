"""テーブルEOSと線形EOSの切り替え(pyR2D2.util.eos_switch)を確認するテスト。

R2D2 本体は相対振幅 ct = max(|ro1|/ro0, |se1|/se0) が ct0 = 3e-3 を
下回るあいだ線形EOSを使い、超えたときだけ表を引く
(runge_kutta_*.F90, artdif_*.F90, cfl.F90, remap_calc.F90 に同じ形で現れる)。
pyR2D2 の従来の eos / Data.eos.eval は常に表を引くので、deep 成層のように
ct が ct0 に届かない計算では出力と食い違う。
"""

import numpy as np
import pytest

from pyR2D2.util import eos_switch, eos_table


CT0 = 3.0e-3


class _Params:
    pass


class _FakeData:
    """eos_table / eos_switch が必要とする属性だけを持つ最小の代役。

    本物の pyR2D2.Data は __getattr__ で self.p へ委譲するので、
    ここでも同じ属性が両方から見えるようにしておく。
    """

    def __init__(self, p):
        self.p = p

    def __getattr__(self, name):
        return getattr(self.p, name)


@pytest.fixture
def data():
    """対数を取ったときに双線形で厳密に表せる、小さな合成EOS表を作る。

    log(q) を (log rho, s) の双線形関数にしておくと、内挿は解析式と
    厳密に一致するので、切り替えの側だけを試験できる。
    """
    nro, nse = 8, 6
    log_ro_e = np.linspace(-20.0, -4.0, nro)
    se_e = np.linspace(1.0e9, 2.0e9, nse)

    p = _Params()
    p.log_ro_e = log_ro_e
    p.se_e = se_e
    p.dlogro_e = log_ro_e[1] - log_ro_e[0]
    p.dse_e = se_e[1] - se_e[0]

    r = log_ro_e[:, None]
    s = (se_e[None, :] - se_e[0]) / (se_e[-1] - se_e[0])
    for var, (a, b, c) in {
        "pr": (10.0, 0.5, 2.0),
        "te": (8.0, 0.2, 1.5),
        "en": (30.0, 0.3, 1.0),
        "op": (-1.0, 0.1, 0.7),
        "dprdro": (25.0, 0.4, 1.2),
    }.items():
        p.__dict__["log_" + var + "_e"] = a + b * r + c * s

    # 鉛直方向の背景。ix = 4 とし、EOS表の範囲内に収める。
    ix = 4
    p.ro0 = np.array([1.0e-6, 1.0e-7, 1.0e-8, 1.0e-9])
    p.se0 = np.linspace(1.2e9, 1.6e9, ix)
    for var in ("pr", "te", "en", "dprdro"):
        p.__dict__[var + "0"] = eos_table(_FakeData(p), p.ro0, p.se0, var)
    # 線形EOSの係数。値そのものに意味は無いので、桁だけ合わせた定数にする。
    p.dprdro, p.dprdse = np.full(ix, 3.0e12), np.full(ix, -1.0e-4)
    p.dtedro, p.dtedse = np.full(ix, 5.0e10), np.full(ix, 2.0e-6)
    p.dendro, p.dendse = np.full(ix, 7.0e12), np.full(ix, 4.0e-4)

    return _FakeData(p)


def _perturbations(data, amplitude):
    """全格子点で ct = amplitude ちょうどになる摂動を作る。"""
    ro0 = data.p.ro0[:, None, None]
    se0 = data.p.se0[:, None, None]
    shape = (data.p.ro0.size, 3, 2)
    ro1 = np.full(shape, amplitude) * ro0
    se1 = np.zeros(shape)
    return ro1, se1


def test_eos_table_matches_scalar_eos(data):
    """配列版 eos_table が、既存のスカラー版と同じ内挿を行う。"""
    from pyR2D2.util import eos

    for var in ("pr", "te", "en", "op", "dprdro"):
        for ro, se in ((1.0e-6, 1.3e9), (3.0e-8, 1.7e9), (2.0e-5, 1.05e9)):
            assert eos_table(data, ro, se, var) == pytest.approx(
                eos(data, ro, se, var), rel=1e-12
            )


def test_below_threshold_uses_linear_eos(data):
    """ct < ct0 では線形EOSそのものになる(表は一切効かない)。"""
    ro1, se1 = _perturbations(data, 0.1 * CT0)
    for var, dro, dse in (
        ("pr", "dprdro", "dprdse"),
        ("te", "dtedro", "dtedse"),
        ("en", "dendro", "dendse"),
    ):
        expected = (
            data.p.__dict__[dro][:, None, None] * ro1
            + data.p.__dict__[dse][:, None, None] * se1
        )
        assert eos_switch(data, ro1, se1, var) == pytest.approx(expected, rel=1e-12)


def test_above_threshold_uses_table(data):
    """ct > ct0 では表の値の摂動になる。"""
    ro1, se1 = _perturbations(data, 10.0 * CT0)
    ro0 = data.p.ro0[:, None, None]
    se0 = data.p.se0[:, None, None]
    for var in ("pr", "te", "en"):
        expected = eos_table(data, ro1 + ro0, se1 + se0, var) - data.p.__dict__[
            var + "0"
        ][:, None, None]
        assert eos_switch(data, ro1, se1, var) == pytest.approx(expected, rel=1e-12)


def test_threshold_is_a_hard_switch(data):
    """閾値をまたぐと不連続に飛ぶ(Fortranの sign() と同じ硬い切り替え)。"""
    below = eos_switch(data, *_perturbations(data, 0.999 * CT0), "pr")
    above = eos_switch(data, *_perturbations(data, 1.001 * CT0), "pr")
    # 摂動そのものは 0.2% しか変わらないのに、値は桁で変わる。
    assert np.abs(above - below).max() > 0.5 * np.abs(below).max()


def test_opacity_is_never_switched(data):
    """op には線形化した対応物が無いので、常に表の値をそのまま返す。"""
    ro0 = data.p.ro0[:, None, None]
    se0 = data.p.se0[:, None, None]
    for amplitude in (0.1 * CT0, 10.0 * CT0):
        ro1, se1 = _perturbations(data, amplitude)
        assert eos_switch(data, ro1, se1, "op") == pytest.approx(
            eos_table(data, ro1 + ro0, se1 + se0, "op"), rel=1e-12
        )


def test_dprdro_falls_back_to_the_background_profile(data):
    """dprdro の線形側は背景プロファイルそのものである(cfl.F90)。"""
    ro1, se1 = _perturbations(data, 0.1 * CT0)
    expected = np.broadcast_to(data.p.dprdro[:, None, None], ro1.shape)
    assert eos_switch(data, ro1, se1, "dprdro") == pytest.approx(expected, rel=1e-12)


def test_vertical_axis_must_come_first(data):
    """第1軸が鉛直でなければ、黙って broadcast せず例外を出す。"""
    bad = np.zeros((data.p.ro0.size + 1, 3, 2))
    with pytest.raises(ValueError):
        eos_switch(data, bad, bad, "pr")


def test_shapes_must_agree(data):
    """ro1 と se1 の形が違えば例外を出す。"""
    ro1, se1 = _perturbations(data, CT0)
    with pytest.raises(ValueError):
        eos_switch(data, ro1, se1[..., :1], "pr")
