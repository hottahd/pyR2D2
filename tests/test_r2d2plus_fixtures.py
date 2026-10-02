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
        np.testing.assert_array_equal(got[k], v, err_msg=k)
    p = Parameters.__new__(Parameters)
    p.read_eos_table(path)
    assert p.ix_e == 7 and p.jx_e == 5
    np.testing.assert_array_equal(p.log_ro_e, want["ro"])
