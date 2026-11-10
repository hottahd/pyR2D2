"""出力の無い番号の time_read は、そのランにある番号を言う (imai、2026-11-10)。

枝分かれしたラン (d008) は主出力が 151 から、τ 面が 1 から始まり、``time_read(150)`` や
``time_read(0, tau=True)`` が ``FileNotFoundError: …/data/time.zip`` で止まっていた。
"""

import shutil
import tarfile
import zipfile
from pathlib import Path

import pytest

import pyR2D2
from pyR2D2.data import _format_steps

FIXTURE = sorted((Path(__file__).parent / "data" / "r2d2plus_fixtures").glob("r2d2plus-*.tar.xz"))[-1]


@pytest.fixture
def branched(tmp_path):
    """見本の legacy ランから 0 番の時刻を消して、1 番から始まる枝分かれのランにする。"""
    with tarfile.open(FIXTURE) as tar:
        tar.extractall(tmp_path, filter="data")
    datadir = tmp_path / "legacy" / "data"
    for family in ("mhd", "tau"):
        (datadir / "time" / family / "t.dac.00000000").unlink()
    return datadir


def test_number_before_the_branch(branched):
    d = pyR2D2.Data(str(branched), verbose=False)
    assert not (branched / "time.zip").exists()
    with pytest.raises(FileNotFoundError, match=r"Output 0 \(mhd\) is not in this run: it has mhd outputs 1-2 "):
        d.time_read(0, verbose=False)
    with pytest.raises(FileNotFoundError, match=r"Output 0 \(tau\).*outputs 1-2.*starts at 1"):
        d.time_read(0, tau=True, verbose=False)
    with pytest.raises(FileNotFoundError, match=r"Output 7 \(mhd\).*outputs 1-2") as err:
        d.time_read(7, verbose=False)
    assert "branched" not in str(err.value)
    assert d.time_read(1, verbose=False) == pyR2D2.Data(str(branched), verbose=False).time_read(1, verbose=False)


def test_number_not_in_time_zip(branched):
    """time.zip があっても、その番号が無ければ KeyError でなく同じ言い方にする。"""
    d = pyR2D2.Data(str(branched), verbose=False)
    t1 = d.time_read(1, verbose=False)
    d.zip_time()
    shutil.rmtree(branched / "time" / "mhd")
    with zipfile.ZipFile(branched / "time.zip") as zf:
        assert "time/mhd/t.dac.00000001" in zf.namelist()
    assert d.time_read(1, verbose=False) == t1
    with pytest.raises(FileNotFoundError, match=r"Output 0 \(mhd\).*outputs 1-2"):
        d.time_read(0, verbose=False)
    with pytest.raises(FileNotFoundError, match=r"Output 0 \(mhd\).*outputs 1-2"):
        d.time_read(0, verbose=False, use_zip=True)


def test_format_steps():
    assert _format_steps(list(range(151, 511))) == "151-510"
    assert _format_steps([1, 2, 3, 7, 9, 10]) == "1-3, 7, 9-10"
    assert _format_steps([1, 3, 5, 7, 9, 11, 13, 15]) == "1, 3, 5, ..., 11, 13, 15"
