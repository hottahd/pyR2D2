"""R2D2plus の各版が書いた出力の見本を読む互換試験 (R2D2plus DEC-606)。

見本 ``tests/data/r2d2plus_fixtures/r2d2plus-<版>.tar.xz`` は、R2D2plus の ``tests/pyr2d2_roundtrip.py --fixture``
が小さな Sod を ``[output] format = "compressed"`` と ``"legacy"`` で走らせて作る。中の ``expected.json`` は、見本を作った
ときの pyR2D2 が読んだ値の SHA-256。pyR2D2 を変えて古い版の出力の読み方が変われば、ここで落ちる。
R2D2plus のタグを打つたびに見本を 1 つ足す (消さない)。
"""

import hashlib
import json
import tarfile
from pathlib import Path

import numpy as np
import pytest

import pyR2D2
from pyR2D2.data_io.restart import RESTART_VARIABLES

FIXTURES = sorted((Path(__file__).parent / "data" / "r2d2plus_fixtures").glob("r2d2plus-*.tar.xz"))
QQ = ["ro", "vx", "vy", "vz", "bx", "by", "bz", "se", "pr", "te", "op"]


def digest(x):
    """R2D2plus tests/pyr2d2_roundtrip.py の array_digest と同じ。"""
    x = np.ascontiguousarray(np.asarray(x))
    h = hashlib.sha256(f"{x.dtype.str}{x.shape}".encode())
    h.update(x.tobytes())
    return h.hexdigest()


@pytest.fixture(params=FIXTURES, ids=[p.name.removesuffix(".tar.xz") for p in FIXTURES])
def fixture_dir(request, tmp_path):
    with tarfile.open(request.param) as tar:
        tar.extractall(tmp_path, filter="data")
    return tmp_path


def test_fixtures_exist():
    assert FIXTURES, "見本が 1 つも無い"


@pytest.mark.parametrize("fmt", ["legacy", "compressed"])
def test_outputs(fixture_dir, fmt):
    if fmt == "compressed":
        pytest.importorskip("zstandard")
    exp = json.loads((fixture_dir / "expected.json").read_text())
    d = pyR2D2.Data(str(fixture_dir / fmt / "data"), verbose=False)
    assert d.p.output_format == fmt
    for n in range(exp["outputs"]):
        d.qf.read(n)
        got = {k: digest(getattr(d.qf, k)) for k in QQ}
        bad = [k for k in QQ if got[k] != exp["qq"][str(n)][k]]
        assert not bad, f"remap/qq {n}: {bad}"
        d.vc.read(n)
        for k, h in exp["vl"][str(n)].items():
            assert digest(getattr(d.vc, k)) == h, f"remap/vl {n} {k}"


@pytest.mark.parametrize("fmt", ["legacy", "compressed"])
def test_restart(fixture_dir, fmt):
    pytest.importorskip("zstandard")  # compressed の restart は圧縮 (legacy のランは無圧縮)
    exp = json.loads((fixture_dir / "expected.json").read_text())["restart"]
    r = pyR2D2.read_restart(fixture_dir / fmt / "data")
    assert r.step == exp["step"] and list(r.shape) == exp["shape"]
    for k in RESTART_VARIABLES:
        assert digest(r.state[k]) == exp["state"][k], k


TAU = ["ro", "se", "pr", "te", "vx", "vy", "vz", "bx", "by", "bz", "fr"]
SLICE = ["ro", "vx", "vy", "vz", "bx", "by", "bz", "se"]


def test_tau_slice_formats_agree(fixture_dir):
    """tau 面と slice を compressed と legacy で読み比べる (tau 面のある見本だけ。v1.1.0-1800 以降)。

    R2D2plus の ctest pyr2d2_roundtrip と同じ規則。高さ he は、compressed が単精度の r - rstar に倍精度で rstar を
    足して返すので単精度の 1 ulp まで違ってよく、面が無い所は legacy が 0・compressed の height が 0。slice の pr・te は
    compressed が単精度の ro・se から計算するので相対 1e-6 まで。2026-10-27 に、astana の 2.4 GB の実ラン
    (io_fixture_20260927) を読んでいた試験をこれに置き換えた。
    """
    pytest.importorskip("zstandard")
    if not any((fixture_dir / "legacy" / "data" / "tau").glob("qq.dac.*")):
        pytest.skip("この見本には tau 面が無い")
    dc = pyR2D2.Data(str(fixture_dir / "compressed" / "data"), verbose=False)
    dl = pyR2D2.Data(str(fixture_dir / "legacy" / "data"), verbose=False)
    found = 0
    for n in range(int(dl.p.nd_tau) + 1):
        dc.qt.read(n)
        dl.qt.read(n)
        for t in ["", "01", "001"]:
            for k in TAU + (["rt"] if not t else []):
                np.testing.assert_array_equal(getattr(dc.qt, k + t), getattr(dl.qt, k + t), err_msg=f"tau {n} {k + t}")
            he, y, h = getattr(dc.qt, "he" + t), getattr(dl.qt, "he" + t), getattr(dc.qt, "height" + t)
            assert np.where(y == 0, h == 0, np.abs(he - y) <= np.spacing(y)).all(), f"tau {n} he{t}"
            found += int(np.count_nonzero(y))
        for direc in ["x", "y", "z"]:
            dc.qs.read(0, direc, n)
            dl.qs.read(0, direc, n)
            for k in SLICE:
                np.testing.assert_array_equal(getattr(dc.qs, k), getattr(dl.qs, k), err_msg=f"slice {direc} {n} {k}")
            for k in ["pr", "te"]:
                y = getattr(dl.qs, k)
                np.testing.assert_allclose(getattr(dc.qs, k), y, rtol=0, atol=1e-6 * np.abs(y).max(),
                                           err_msg=f"slice {direc} {n} {k}")
    assert found > 0, "tau 面がどこにも見つかっていない (比べたのが 0 だけ)"


def test_checksums_and_initial_state_from_output(fixture_dir, tmp_path):
    """R2D2plus が書いた .z の XXH3-64 を照合でき、legacy の出力を初期条件 (write_initial_state) に書き直せる。

    2026-10-27 まで astana の io_fixture_20260927 (128^3、4 ランク) で行っていた試験を見本に移した。
    """
    pytest.importorskip("xxhash")
    from pyR2D2.data_io.compressed import read_r2d2plus_z

    zfiles = sorted((fixture_dir / "compressed" / "data" / "remap" / "qq").rglob("qq.z.*"))
    assert zfiles
    for f in zfiles:
        read_r2d2plus_z(f, verify=True)  # xxhash の XXH3-64 が R2D2plus (XXH3_64bits) と同じ
    d = pyR2D2.Data(str(fixture_dir / "legacy" / "data"), verbose=False)
    keys = ["ro", "vx", "vy", "vz", "bx", "by", "bz", "se"]
    d.qf.read(0, keys=keys)
    fields = {k: getattr(d.qf, k) for k in keys}
    path = pyR2D2.write_initial_state(tmp_path / "init.z", fields, attributes={"nd": 0})
    _, arrays = read_r2d2plus_z(path, verify=True)
    for k in keys:
        assert arrays[k].shape == (d.ix, d.jx, d.kx)
        np.testing.assert_array_equal(arrays[k], getattr(d.qf, k).astype(np.float64), err_msg=k)


# R2D2plus-input の書き手 (contracts/rtable.write_table) が書いた r2d2plus-table。
# R2D2plus-input の tests/make_pyr2d2_table_fixture.py が作る。配列はそちらと同じ式で作り直して比べる。
TABLES = sorted((Path(__file__).parent / "data" / "r2d2plus_fixtures").glob("rtable-*.tbl"))
TABLE_NAMES = ["pr", "en", "te", "op", "dprdro", "cv", "cp", "dprdse", "dtedro", "dtedse", "dendro", "dendse"]


def table_arrays():
    ro = np.linspace(-12.0, -2.0, 7)
    se = np.linspace(1.0e9, 7.0e9, 5)
    out = {"ro": ro, "se": se}
    for i, name in enumerate(TABLE_NAMES):
        out[name] = np.exp(np.add.outer((0.1 + 0.01 * i) * ro, (1.0 + 0.1 * i) * 1.0e-10 * se))
    return out


@pytest.mark.parametrize("path", TABLES, ids=[p.name for p in TABLES])
def test_rtable_written_by_r2d2plus_input(path):
    from pyR2D2.data_io.compressed import load_eos_table_raw
    from pyR2D2.data_io.parameters import Parameters

    want = table_arrays()
    got = load_eos_table_raw(path)
    assert set(got) == set(want)
    for k, v in want.items():
        # np.exp は機械で 1 ulp 違う (GitHub の CI で 1.4e-16)。読み違いなら桁ごと違うので相対 1e-14 で比べる
        np.testing.assert_allclose(got[k], v, rtol=1e-14, atol=0, err_msg=k)
    p = Parameters.__new__(Parameters)
    p.read_eos_table(path)
    assert p.ix_e == 7 and p.jx_e == 5
    np.testing.assert_array_equal(p.log_ro_e, want["ro"])
