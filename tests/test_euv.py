"""R2D2plus の量の並び 2 (DEC-671、2026-10-31) の読み手と、EM から強度を作る道具 (pyR2D2.synth.euv) の試験。

- vl_hist.dac・vl_dem.dac を OnTheFly.read_hist・read_dem が正しい形と区画で読むこと、無ければ None
- EOS 表の任意の配列 ne を read_eos_table が読むこと
- ResponseTable の保存と読み込み、箱の平均、intensity の和
- em_from_state が vl_dem と同じ定義 (表の T と ne、鉛直は格子幅、箱の外は数えない) で EM を作ること
- 同梱の AIA・XRT の表が読めて、チャンネルと形がそろっていること
"""

from types import SimpleNamespace

import numpy as np
import pytest

from pyR2D2.data_io.parameters import Parameters
from pyR2D2.data_io.read import OnTheFly
from pyR2D2.synth import euv


def _bare_onthefly(tmp_path, ix, jx, kx):
    vc = OnTheFly.__new__(OnTheFly)
    vc.data = None
    vc.datadir = tmp_path
    vc.endian = "<"
    vc.ix, vc.jx, vc.kx = ix, jx, kx
    return vc


def test_read_hist_and_dem(tmp_path):
    ix, jx, kx, n = 6, 4, 3, 7
    vldir = tmp_path / "remap" / "vl"
    vldir.mkdir(parents=True)
    hist = np.arange(ix * 180, dtype="<f4").reshape((ix, 180), order="F")
    hist.flatten(order="F").tofile(vldir / f"vl_hist.dac.{n:08d}")
    # 3 区画: z 方向の視線 (ix, jx, 60)、y 方向 (ix, kx, 60)、鉛直 (jx, kx, 60)
    secs = [np.full((ix, jx, 60), 1, "<f4"), np.full((ix, kx, 60), 2, "<f4"), np.full((jx, kx, 60), 3, "<f4")]
    secs[0][1, 2, 3] = 10.0
    secs[2][0, 1, 59] = 30.0
    np.concatenate([s.flatten(order="F") for s in secs]).tofile(vldir / f"vl_dem.dac.{n:08d}")

    vc = _bare_onthefly(tmp_path, ix, jx, kx)
    vc.read_hist(n)
    vc.read_dem(n)
    np.testing.assert_array_equal(vc.hist_te, hist[:, :100])
    np.testing.assert_array_equal(vc.hist_ro, hist[:, 100:])
    assert vc.em_z.shape == (ix, jx, 60) and vc.em_y.shape == (ix, kx, 60) and vc.em_x.shape == (jx, kx, 60)
    np.testing.assert_array_equal(vc.em_z, secs[0])
    np.testing.assert_array_equal(vc.em_y, secs[1])
    np.testing.assert_array_equal(vc.em_x, secs[2])
    assert vc.hist_logte_edges.size == 101 and vc.hist_logro_edges.size == 81
    np.testing.assert_allclose(vc.dem_logte_edges[[0, -1]], [5.0, 8.0])

    # 無ければ None (前の値は消す)
    vc.read_hist(n + 1)
    vc.read_dem(n + 1)
    assert vc.hist_te is None and vc.hist_ro is None and vc.em_x is None and vc.em_z is None

    # 寸法の合わないファイルは黙って読まない
    np.zeros(10, "<f4").tofile(vldir / f"vl_dem.dac.{n + 2:08d}")
    with pytest.raises(ValueError, match="expected"):
        vc.read_dem(n + 2)


def _eos_npz(path, with_ne):
    log_ro = np.linspace(-30.0, 0.0, 31)
    se = np.linspace(0.0, 3.0e9, 31)
    one = np.ones((31, 31))
    d = dict(ro=log_ro, se=se, pr=one, en=one, te=one, op=one, dprdro=one)
    if with_ne:
        d["ne"] = 2.0 * one
    np.savez(path, **d)


@pytest.mark.parametrize("with_ne", [True, False])
def test_read_eos_table_ne(tmp_path, with_ne):
    path = tmp_path / "eos.npz"
    _eos_npz(path, with_ne)
    p = Parameters.__new__(Parameters)
    p.log_ne_e = "old"
    p.read_eos_table(path)
    if with_ne:
        np.testing.assert_allclose(p.log_ne_e, np.log(2.0))
    else:
        assert "log_ne_e" not in p.__dict__


def test_response_table_roundtrip_and_intensity(tmp_path):
    logt = np.linspace(4.0, 8.0, 81)
    table = euv.ResponseTable(logt, {"A": np.ones_like(logt), "B": logt.copy()}, meta={"source": "test"})
    for name in ("t.txt", "t.npz"):
        back = euv.ResponseTable.load(table.save(tmp_path / name) or tmp_path / name)
        assert back.channels == ["A", "B"] and back.meta["source"] == "test"
        np.testing.assert_allclose(back.logt, logt)
        np.testing.assert_allclose(back["B"], logt, rtol=1e-6)
    # 箱の平均: 定数はそのまま、log T に線形なら箱の中心の値
    edges = euv.DEM_LOGTE_EDGES
    np.testing.assert_allclose(table.bin_mean("A", edges), 1.0)
    np.testing.assert_allclose(table.bin_mean("B", edges), 0.5 * (edges[1:] + edges[:-1]))
    # 表の外は 0
    out = euv.ResponseTable(np.array([6.0, 7.0]), {"A": np.array([1.0, 1.0])})
    m = out.bin_mean("A", edges)
    assert m[0] == 0.0 and m[-1] == 0.0 and m[30] == 1.0

    em = np.zeros((2, 3, 60))
    em[0, 1, 10] = 3.0
    em[1, 2, :] = 1.0
    img = euv.intensity(em, table)
    assert img["A"].shape == (2, 3)
    assert img["A"][0, 1] == pytest.approx(3.0)
    assert img["A"][1, 2] == pytest.approx(60.0)
    assert img["B"][0, 1] == pytest.approx(3.0 * (5.0 + 0.05 * 10.5))
    with pytest.raises(ValueError, match="bins"):
        euv.intensity(np.zeros((2, 59)), table)
    # ファイル名でも渡せる
    table.save(tmp_path / "x.txt")
    assert euv.intensity(em, tmp_path / "x.txt", channels=["A"]).keys() == {"A"}


class _FakeQF:
    def __init__(self, ro, se):
        self._ro, self._se = ro, se

    def read(self, n, keys="all"):
        self.ro, self.se = self._ro, self._se


def _fake_data(ix=8, jx=3, kx=4, seed=1):
    """log T = a log ρ + b s、log ne = c log ρ の表 (双線形の補間で厳密) を持つ偽のラン。"""
    rng = np.random.default_rng(seed)
    log_ro_e = np.linspace(-40.0, 0.0, 41)
    se_e = np.linspace(0.0, 1.0, 11)
    log10_te = 4.5 + np.add.outer(-0.08 * (log_ro_e + 30.0), 0.0 * se_e) + 2.0 * se_e[None, :]
    log_ne = np.add.outer(0.9 * log_ro_e + 60.0, 0.0 * se_e)
    p = SimpleNamespace(
        log_ro_e=log_ro_e, se_e=se_e, dlogro_e=1.0, dse_e=0.1,
        log_te_e=np.log(10.0) * log10_te, log_ne_e=log_ne,
        ro0=np.exp(np.linspace(-35.0, -25.0, ix)), se0=np.full(ix, 0.3),
        geometry="Cartesian", margin=2,
        # 鉛直は不均等 (両端の幅が片側差分と違うように)、水平は等間隔
        xg=np.cumsum(np.linspace(1.0e6, 3.0e6, ix + 4) ** 1.3),
        yg=np.arange(jx + 4) * 2.0e6, zg=np.arange(kx + 4) * 5.0e6,
    )
    ro = (p.ro0[:, None, None] * rng.uniform(-0.5, 0.5, (ix, jx, kx))).astype("f4")
    se = rng.uniform(-0.3, 0.6, (ix, jx, kx)).astype("f4")
    data = SimpleNamespace(p=p, qf=_FakeQF(ro, se), log_ro_e=log_ro_e, se_e=se_e, dlogro_e=1.0, dse_e=0.1)
    return data


def test_em_from_state_matches_definition():
    data = _fake_data()
    p = data.p
    ro = data.qf._ro.astype(np.float64) + p.ro0[:, None, None]
    se = data.qf._se.astype(np.float64) + p.se0[:, None, None]
    lr = np.log(ro)
    log10_te = 4.5 - 0.08 * (lr + 30.0) + 2.0 * se
    ne = np.exp(0.9 * lr + 60.0)
    b = np.floor((log10_te - 5.0) / 0.05).astype(int)
    ix, jx, kx = ro.shape
    xg = p.xg
    dx = np.array([0.5 * (xg[i + 3] - xg[i + 1]) for i in range(ix)])
    assert dx[0] != xg[3] - xg[2]  # 端も margin を使った中心差分 (np.gradient の片側差分とは違う)
    ref = {"x": np.zeros((jx, kx, 60)), "y": np.zeros((ix, kx, 60)), "z": np.zeros((ix, jx, 60))}
    inside = 0
    for i in range(ix):
        for j in range(jx):
            for k in range(kx):
                if 0 <= b[i, j, k] < 60:
                    inside += 1
                    w = ne[i, j, k] ** 2
                    ref["x"][j, k, b[i, j, k]] += w * dx[i]
                    ref["y"][i, k, b[i, j, k]] += w * 2.0e6
                    ref["z"][i, j, b[i, j, k]] += w * 5.0e6
    assert 0 < inside < ix * jx * kx  # 箱の内と外の両方がある
    for los in "xyz":
        np.testing.assert_allclose(euv.em_from_state(data, 0, los), ref[los], rtol=1e-9)
    with pytest.raises(ValueError, match="los"):
        euv.em_from_state(data, 0, "w")
    del p.log_ne_e
    with pytest.raises(ValueError, match="electron density"):
        euv.em_from_state(data, 0, "x")


@pytest.mark.parametrize("name, channels", [("aia", euv.AIA_CHANNELS), ("xrt", euv.XRT_CHANNELS)])
def test_bundled_tables(name, channels):
    table = euv.load_response(name)
    assert tuple(table.channels) == tuple(channels)
    assert table.units == "DN cm^5 s^-1 pix^-1"
    assert table.logt[0] <= 5.0 and table.logt[-1] >= 8.0 - 1e-9
    for ch in channels:
        assert np.all(table[ch] >= 0.0) and table[ch].max() > 0.0
    assert "source" in table.meta and "url" in table.meta


def test_bundled_aia_peaks():
    """AIA の主な山の温度 (Boerner et al. 2012 の図の位置)。"""
    t = euv.load_response("aia")
    peak = {ch: t.logt[np.argmax(np.where(t.logt > 5.5, t[ch], 0))] for ch in t.channels}
    assert 5.8 <= peak["A171"] <= 5.95
    assert 6.1 <= peak["A193"] <= 6.3
    assert 6.2 <= peak["A211"] <= 6.4
    assert 6.7 <= peak["A94"] <= 7.0
