"""出力の間隔が 4 つのラン (R2D2plus output_cadence = 2、2026-10-31) と、three_d_every のランの試験。

R2D2plus の見本 (tests/data/r2d2plus_fixtures) を展開し、数え方 2 のランの形 (nd.dac の 4 列、params.dac の
output_cadence・dtout_vl・dtout_slice、time/vl・time/slice、compressed なら format.toml の版 3 と {nd_slice} の雛形)
に書き換えて読む。
"""

import shutil
import tarfile
import zipfile
from pathlib import Path

import numpy as np
import pytest

import pyR2D2

FIXTURE = sorted((Path(__file__).parent / "data" / "r2d2plus_fixtures").glob("r2d2plus-*.tar.xz"))[-1]


def extract(dest, fmt):
    if fmt == "compressed":
        try:
            import numcodecs  # noqa: F401
        except ImportError:
            pytest.importorskip("zstandard")
    with tarfile.open(FIXTURE) as tar:
        tar.extractall(dest, filter="data")
    return dest / fmt / "data"


def write_time(path, t, endian="<"):
    path.parent.mkdir(parents=True, exist_ok=True)
    np.array([t], dtype=endian + "d").tofile(path)


def to_cadence2(datadir, fmt):
    """見本のランを数え方 2 の形にする。vl・slice の時刻は 3D・tau と違う値にして、取り違えを見分ける。"""
    param = datadir / "param"
    (param / "nd.dac").write_text("       2       2       2       2\n")
    with open(param / "params.dac", "a") as f:
        f.write("         2 output_cadence i\n")
        f.write("0.2500000000E+00 dtout_vl d\n")
        f.write("0.2500000000E+00 dtout_slice d\n")
    for n in range(3):
        write_time(datadir / "time" / "vl" / f"t.dac.{n:08d}", 100.0 + n)
        write_time(datadir / "time" / "slice" / f"t.dac.{n:08d}", 200.0 + n)
    if fmt == "compressed":
        text = (param / "format.toml").read_text()
        text = text.replace("version = 1", "version = 3").replace("{nd_tau:08d}.{n:08d}", "{nd_slice:08d}.{n:08d}")
        text += "output_cadence = 2\ndtout = 0.5\ndtout_vl = 0.25\ndtout_tau = 0.5\ndtout_slice = 0.25\n"
        (param / "format.toml").write_text(text)


@pytest.mark.parametrize("fmt", ["legacy", "compressed"])
def test_old_cadence_maps_vl_and_slice(tmp_path, fmt):
    d = pyR2D2.Data(str(extract(tmp_path, fmt)), verbose=False)
    assert d.p.output_cadence == 1
    assert d.p.nd_vl == d.p.nd and d.p.nd_slice == d.p.nd_tau
    assert d.p.dtout_vl == d.p.dtout and d.p.dtout_slice == d.p.dtout_tau
    for n in range(3):
        assert d.time_read(n, kind="vl", verbose=False) == d.time_read(n, verbose=False)
        assert d.time_read(n, kind="slice", verbose=False) == d.time_read(n, tau=True, verbose=False)
    assert d.time_steps("vl") == [0, 1, 2]
    assert d.qq_steps() == [0, 1, 2]


@pytest.mark.parametrize("fmt", ["legacy", "compressed"])
def test_cadence2_reads_own_times(tmp_path, fmt, capsys):
    datadir = extract(tmp_path, fmt)
    orig = pyR2D2.Data(str(datadir), verbose=False)
    orig.qs.read(1, "x", 2)
    expected = orig.qs.ro.copy()
    to_cadence2(datadir, fmt)
    d = pyR2D2.Data(str(datadir), verbose=False)
    assert d.p.output_cadence == 2
    assert (d.p.nd, d.p.nd_vl, d.p.nd_tau, d.p.nd_slice) == (2, 2, 2, 2)
    assert d.p.dtout_vl == 0.25 and d.p.dtout_slice == 0.25
    assert d.time_read(1, kind="vl", verbose=False) == 101.0
    assert d.time_read(2, kind="slice", verbose=False) == 202.0
    assert d.time_read(1, verbose=False) == orig.time_read(1, verbose=False)
    steps, t = d.times("slice")
    assert list(steps) == [0, 1, 2] and list(t) == [200.0, 201.0, 202.0]
    # slice のファイルの番号は nd_slice (compressed の版 3 は雛形も {nd_slice})
    d.qs.read(1, "x", 2)
    np.testing.assert_array_equal(d.qs.ro, expected)
    d.vc.read(2)
    d.p.summary()
    out = capsys.readouterr().out
    assert "vl    nd_vl = 2  dtout_vl = 0.25 s" in out
    assert "slice nd_slice = 2  dtout_slice = 0.25 s" in out


def test_nd_dac_minus_one_and_time_dirs(tmp_path):
    datadir = extract(tmp_path, "legacy")
    to_cadence2(datadir, "legacy")
    (datadir / "param" / "nd.dac").write_text("2 2 -1 -1\n")
    d = pyR2D2.Data(str(datadir), verbose=False)
    assert d.p.nd_vl == 2 and d.p.nd_slice == 2  # time/vl・time/slice の番号から
    shutil.rmtree(datadir / "time" / "vl")
    shutil.rmtree(datadir / "time" / "slice")
    d = pyR2D2.Data(str(datadir), verbose=False)
    assert d.p.nd_vl == -1 and d.p.nd_slice == -1


def set_three_d_every(datadir, k, seg=None):
    runs = sorted((datadir / "param" / "runs").iterdir())
    seg = seg or runs[0]
    path = seg / "run.toml"
    text = path.read_text()
    text = text.replace("[output]\n", f"[output]\nthree_d_every = {k}\n", 1)
    path.write_text(text)


def test_qq_steps_three_d_every(tmp_path, capsys):
    datadir = extract(tmp_path, "legacy")
    for n in range(9):
        write_time(datadir / "time" / "mhd" / f"t.dac.{n:08d}", 0.5 * n)
    (datadir / "param" / "nd.dac").write_text("8 2\n")
    set_three_d_every(datadir, 2)
    d = pyR2D2.Data(str(datadir), verbose=False)
    assert d.three_d_every_values() == [2]
    assert d.qq_steps() == [0, 1, 2, 4, 6, 8]  # 0, 2, 4, 6, 8 と、ファイルのある 1
    d.p.summary()
    assert "3D (qq) steps = 6 (three_d_every = 2)" in capsys.readouterr().out

    # 2 つめの区間 (時刻 1.6 から再開) は three_d_every = 3: nd = 4..8 のうち 6
    first = sorted((datadir / "param" / "runs").iterdir())[0]
    second = first.parent / "20991231-000000"
    shutil.copytree(first, second)
    text = (second / "run.toml").read_text().replace("three_d_every = 2", "three_d_every = 3")
    (second / "run.toml").write_text(text)
    prov = second / "provenance.toml"
    lines = [l for l in prov.read_text().splitlines() if not l.startswith("resume_time")]
    i = lines.index("[segment]")
    lines.insert(i + 1, 'resume_time = "1.6"')
    prov.write_text("\n".join(lines) + "\n")
    d = pyR2D2.Data(str(datadir), verbose=False)
    assert d.three_d_every_values() == [2, 3]
    assert d.qq_steps() == [0, 1, 2, 6]
    # restart の時刻が出力の時刻ちょうど (nd = 3 の 1.5) なら、nd = 3 は前の区間のもの。後の区間は 4 から
    prov.write_text(prov.read_text().replace('resume_time = "1.6"', 'resume_time = "1.5000000000000000e+00"'))
    assert pyR2D2.Data(str(datadir), verbose=False).qq_steps() == [0, 1, 2, 6]
    # 後の区間が 3 から始まるなら (時刻 1.4 から再開) 3・6 が 3D
    prov.write_text(prov.read_text().replace('resume_time = "1.5000000000000000e+00"', 'resume_time = "1.4"'))
    assert pyR2D2.Data(str(datadir), verbose=False).qq_steps() == [0, 1, 2, 3, 6]


def test_time_zip_member_names(tmp_path):
    """time.zip の中は time/<族>/… (zip_time)。2026-10-31 までの r2d2plus-archive の <族>/… も読める。"""
    datadir = extract(tmp_path, "legacy")
    t1 = pyR2D2.Data(str(datadir), verbose=False).time_read(1, verbose=False)
    with zipfile.ZipFile(datadir / "time.zip", "w") as zf:
        for f in (datadir / "time").rglob("t.dac.*"):
            zf.write(f, str(f.relative_to(datadir / "time")))
    shutil.rmtree(datadir / "time")
    d = pyR2D2.Data(str(datadir), verbose=False)
    assert d.time_read(1, verbose=False) == t1
    assert d.time_steps("mhd") == [0, 1, 2]
