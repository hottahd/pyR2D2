"""R2D2plus (C++版) 由来の追随の読み込みテスト。

- vl_et.dac (C++独自の ⟨ρe⟩、2026-08-31追加) を OnTheFly.read_et が読めること、
  無い場合 (Fortranラン) に etm = None で黙って戻ること。
- Parameters.read_eos_table が明示パスで任意寸法の npz (A2拡張表 177x166 など)
  を読めること、省略時にファイルが無ければ従来どおり何もしないこと。
"""

import numpy as np

from pyR2D2.data_io.parameters import Parameters
from pyR2D2.data_io.read import OnTheFly


def _bare_onthefly(tmp_path, ix, jx):
    """__init__ を通さず read_et に必要な属性だけ持つ OnTheFly を作る。"""
    vc = OnTheFly.__new__(OnTheFly)
    vc.data = None  # __getattr__ の委譲先 (触られない)
    vc.datadir = tmp_path
    vc.endian = "<"
    vc.ix = ix
    vc.jx = jx
    return vc


def test_read_et_reads_r2d2plus_file(tmp_path):
    ix, jx, n = 5, 3, 12
    vldir = tmp_path / "remap" / "vl"
    vldir.mkdir(parents=True)
    # 書き込み側 (legacy_on_the_fly_writer.cpp の write_et) と同じ
    # 単精度・Fortran順 (鉛直 ix が最速)
    expected = np.arange(ix * jx, dtype="<f4").reshape((ix, jx), order="F")
    expected.flatten(order="F").tofile(vldir / f"vl_et.dac.{n:08d}")

    vc = _bare_onthefly(tmp_path, ix, jx)
    vc.read_et(n)

    assert vc.etm is not None
    assert vc.etm.shape == (ix, jx)
    np.testing.assert_array_equal(vc.etm, expected)


def test_read_et_missing_file_sets_none(tmp_path):
    (tmp_path / "remap" / "vl").mkdir(parents=True)
    vc = _bare_onthefly(tmp_path, 4, 4)
    vc.etm = np.zeros((4, 4))  # 前の時刻の値が残っていても上書きされること
    vc.read_et(1)
    assert vc.etm is None


def _synthetic_eos_npz(path, nro, nse):
    log_ro = np.linspace(-41.0, 0.0, nro)
    se = np.linspace(9.0e8, 7.5e9, nse)
    pr = np.exp(np.add.outer(0.1 * log_ro, 1.0e-9 * se))
    np.savez(
        path,
        ro=log_ro,
        se=se,
        pr=pr,
        en=2.0 * pr,
        te=3.0 * pr,
        op=4.0 * pr,
        dprdro=5.0 * pr,
    )
    return log_ro, se, pr


def test_read_eos_table_explicit_path_any_shape(tmp_path):
    # A2拡張表と同じく正方でない寸法 (R2D2plus DEC-279/283 は 177x166)
    nro, nse = 7, 5
    npz = tmp_path / "eos_table_sero_a2ext.npz"
    log_ro, se, pr = _synthetic_eos_npz(npz, nro, nse)

    p = Parameters.__new__(Parameters)
    p.read_eos_table(npz)

    assert p.ix_e == nro
    assert p.jx_e == nse
    np.testing.assert_array_equal(p.log_ro_e, log_ro)
    np.testing.assert_allclose(p.log_pr_e, np.log(pr + 1.0e-200))
    np.testing.assert_allclose(p.dlogro_e, log_ro[1] - log_ro[0])
    np.testing.assert_allclose(p.dse_e, se[1] - se[0])


def test_read_eos_table_default_missing_is_silent(tmp_path):
    # C++ランの形: <run>/data はあるが <run>/input_data が無い
    datadir = tmp_path / "data"
    datadir.mkdir()
    p = Parameters.__new__(Parameters)
    p.datadir = datadir
    p.read_eos_table()
    assert not hasattr(p, "log_ro_e")
