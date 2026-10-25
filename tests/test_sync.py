"""pyR2D2.Sync の試験。

リモートの ``~/work/<project>/run/<caseid>/`` を手元の一時ディレクトリに作り、rsync の
``server:path`` を手元のパスに読み替えて実際に rsync を走らせる。従来形式と R2D2plus の
新出力形式 (DEC-576)・restart (DEC-076) で、欲しいファイルだけが来ること、
``data/param/source`` と ``data/param/bin`` (DEC-638) が来ないことを確かめる。
"""

import glob
import shutil
import subprocess
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from pyR2D2.data_io.compressed import _FORMAT_DEFAULTS
from pyR2D2.sync.sync import Sync

pytestmark = pytest.mark.skipif(shutil.which("rsync") is None, reason="needs rsync")

PROJECT, CASEID, SERVER = "proj", "d001", "remote"
NPE = 2


def touch(path, text="x"):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)


@pytest.fixture
def remote(tmp_path, monkeypatch):
    """リモートのホーム (``tmp_path/home``) と、そこを読む rsync / ssh の差し替え。"""
    home = tmp_path / "home"
    home.mkdir()

    def to_local(arg):
        if arg.startswith(SERVER + ":"):
            path = arg[len(SERVER) + 1 :]
            # rsync では末尾の / に意味がある (中身を写す) ので落とさない
            trailing = "/" if path.endswith("/") else ""
            path = path if path.startswith("/") else str(home / path)
            # リモートのシェルがするワイルドカードの展開
            if any(c in path for c in "*?["):
                return [p + trailing for p in sorted(glob.glob(path))]
            return [path.rstrip("/") + trailing]
        return [arg]

    def fake_rsync(args):
        local = []
        skip = False
        for a in args:
            if skip:
                skip = False
                continue
            if a == "-e":
                skip = True
                continue
            local += to_local(a)
        return subprocess.run(["rsync", "-a"] + local, capture_output=True, text=True)

    def fake_exists(server, path, ssh="ssh"):
        return bool(glob.glob(str(home / path)))

    monkeypatch.setattr(Sync, "rsync_subprocess_wrapper", staticmethod(fake_rsync))
    monkeypatch.setattr(Sync, "_remote_exists", staticmethod(fake_exists))
    return home / "work" / PROJECT / "run" / CASEID


def make_sync(tmp_path, compressed, nd=3):
    datadir = tmp_path / "local" / PROJECT / "run" / CASEID / "data"
    datadir.mkdir(parents=True)
    format_info = dict(_FORMAT_DEFAULTS, format="compressed") if compressed else None
    p = SimpleNamespace(
        output_format="compressed" if compressed else "legacy",
        format_info=format_info,
        np_ijr=np.array([[0], [1]]),
        x=np.array([0.0, 1.0]),
        i2ir=np.array([1, 2]),
        nd=nd,
        npe=NPE,
    )
    data = SimpleNamespace(datadir=datadir, p=p, **vars(p))
    return Sync(data), datadir


def local_files(datadir):
    return sorted(str(p.relative_to(datadir)) for p in datadir.rglob("*") if p.is_file())


def make_compressed_remote(run):
    data = run / "data"
    touch(data / "param" / "format.toml", 'format = "compressed"\n')
    touch(data / "param" / "params.dac")
    touch(data / "param" / "source" / "0123.r2d2plus_source.tar.zst")
    touch(data / "param" / "bin" / "0123" / "r2d2plus_run")
    touch(data / "param" / "inputs" / "0123_star_def.toml")
    touch(data / "time" / "mhd" / "t.dac.00000001")
    touch(data / "time" / "tau" / "t.dac.00000005")
    for nd in [1, 2]:
        for rank in range(NPE):
            touch(data / f"remap/qq/00000/{rank:08d}/qq.z.{nd:08d}.{rank:08d}")
    touch(data / "remap" / "qq" / "qq.c.00000001")
    touch(data / "remap" / "vl" / "vl_xy.dac.00000001")
    touch(data / "tau" / "qq.c.00000005")
    touch(data / "tau" / "qq.c.00000006")
    touch(data / "slice" / "params.dac")
    touch(data / "slice" / "slice.dac")
    for direc in "xyz":
        for n in [1, 2]:
            touch(data / "slice" / f"qq{direc}.c.00000005.{n:08d}")
    touch(data / "prev" / "00000" / "00000000" / "qq.z.00000001.00000001.00000000")
    touch(data / "photosphere" / "meta.toml")
    # restart: v2 の交代 slot (e が古い・o が新しい) と恒久 slot、v1 の恒久 slot
    for slot, step in [("e", 10), ("o", 20), ("00000002", 15)]:
        touch(
            data / "restart" / slot / "meta.toml",
            f'[restart]\nformat_version = 2\ngeneration = ".generation-new"\nstep = {step}\n',
        )
        for gen in [".generation-new", ".generation-old"]:
            for rank in range(NPE):
                touch(data / "restart" / slot / gen / f"state.{rank:08d}.bin")
            touch(data / "restart" / slot / gen / "extra.afc.bin")
    touch(data / "restart" / "00000001" / "meta.toml", "[restart]\nformat_version = 1\nstep = 5\n")
    for rank in range(NPE):
        touch(data / "restart" / "00000001" / f"state.{rank:08d}.bin")
    touch(run / "run.toml")
    touch(run / "input_data" / "star_def.toml")


def test_setup_skips_bulk_and_source(tmp_path, remote, monkeypatch):
    make_compressed_remote(remote)
    dist = tmp_path / "local" / PROJECT / "run"
    dist.mkdir(parents=True)
    Sync.setup(SERVER, CASEID, project=PROJECT, dist=str(dist) + "/")
    got = local_files(dist / CASEID)
    assert "data/param/format.toml" in got
    assert "data/param/inputs/0123_star_def.toml" in got
    assert "data/slice/params.dac" in got and "data/slice/slice.dac" in got
    assert "run.toml" in got and "input_data/star_def.toml" in got
    for unwanted in ["param/source", "param/bin", "restart", "photosphere", "remap/qq",
                     "tau/qq", "slice/qq", "prev"]:
        assert not any(unwanted in f for f in got), unwanted


def test_all_skips_source(tmp_path, remote):
    make_compressed_remote(remote)
    dist = tmp_path / "local" / PROJECT / "run"
    sync, _ = make_sync(tmp_path, compressed=True)
    sync.all(SERVER, project=PROJECT, dist=str(dist) + "/")
    got = local_files(dist / CASEID)
    assert "data/restart/o/meta.toml" in got
    assert "data/param/inputs/0123_star_def.toml" in got
    assert not any("param/source" in f or "param/bin" in f for f in got)


def test_compressed_remap_qq_and_xselect(tmp_path, remote):
    make_compressed_remote(remote)
    sync, datadir = make_sync(tmp_path, compressed=True)
    sync.remap_qq(SERVER, 2, project=PROJECT)
    assert local_files(datadir) == [
        "remap/qq/00000/00000000/qq.z.00000002.00000000",
        "remap/qq/00000/00000001/qq.z.00000002.00000001",
    ]
    sync.xselect(0.0, SERVER, project=PROJECT)  # np_ijr[0, :] = [0]、全ステップ
    assert "remap/qq/00000/00000000/qq.z.00000001.00000000" in local_files(datadir)
    assert "remap/qq/00000/00000001/qq.z.00000001.00000001" not in local_files(datadir)


def test_compressed_tau(tmp_path, remote):
    make_compressed_remote(remote)
    sync, datadir = make_sync(tmp_path, compressed=True)
    sync.tau(SERVER, 6, project=PROJECT)
    assert local_files(datadir) == ["tau/qq.c.00000006"]
    sync.tau(SERVER, project=PROJECT)
    assert local_files(datadir) == [
        "tau/qq.c.00000005",
        "tau/qq.c.00000006",
        "time/tau/t.dac.00000005",
    ]


def test_compressed_slice(tmp_path, remote):
    make_compressed_remote(remote)
    sync, datadir = make_sync(tmp_path, compressed=True)
    sync.slice(SERVER, n=5, direc="y", n_slice=2, project=PROJECT)
    assert local_files(datadir) == [
        "slice/params.dac",
        "slice/qqy.c.00000005.00000002",
        "slice/slice.dac",
    ]
    sync.slice(SERVER, project=PROJECT)
    assert len([f for f in local_files(datadir) if ".c." in f]) == 6


def test_restart_end_step_takes_newer_slot_and_current_generation(tmp_path, remote):
    make_compressed_remote(remote)
    sync, datadir = make_sync(tmp_path, compressed=True)
    sync.check(SERVER, 0, project=PROJECT, end_step=True)
    got = [f for f in local_files(datadir) if f.startswith("restart/o")]
    assert got == [
        "restart/o/.generation-new/extra.afc.bin",
        "restart/o/.generation-new/state.00000000.bin",
        "restart/o/.generation-new/state.00000001.bin",
        "restart/o/meta.toml",
    ]
    assert not any(".generation-new/state" in f for f in local_files(datadir) if f.startswith("restart/e"))


def test_restart_numbered_slots(tmp_path, remote):
    make_compressed_remote(remote)
    sync, datadir = make_sync(tmp_path, compressed=True)
    sync.check(SERVER, 2, project=PROJECT)
    sync.check(SERVER, 1, project=PROJECT)  # v1
    got = local_files(datadir)
    assert "restart/00000002/.generation-new/state.00000001.bin" in got
    assert not any("00000002/.generation-old" in f for f in got)
    assert "restart/00000001/state.00000001.bin" in got


def test_legacy_paths_unchanged(tmp_path, remote):
    data = remote / "data"
    for rank in range(NPE):
        touch(data / f"remap/qq/00000/{rank:08d}/qq.dac.00000000.{rank:08d}")
        touch(data / f"remap/qq/00000/{rank:08d}/qq.dac.00000003.{rank:08d}")
        touch(data / f"qq/00000/{rank:08d}/qq.dac.o.{rank:08d}")
    touch(data / "tau" / "qq.dac.00000004")
    touch(data / "time" / "tau" / "t.dac.00000004")
    touch(data / "slice" / "params.dac")
    touch(data / "slice" / "slice.dac")
    touch(data / "slice" / "qqz.dac.00000004.00000001")
    sync, datadir = make_sync(tmp_path, compressed=False, nd=3)
    sync.remap_qq(SERVER, 3, project=PROJECT)
    sync.tau(SERVER, project=PROJECT)
    sync.slice(SERVER, n=4, direc="z", n_slice=1, project=PROJECT)
    sync.check(SERVER, 0, project=PROJECT, end_step=True)
    assert local_files(datadir) == [
        "qq/00000/00000000/qq.dac.o.00000000",
        "qq/00000/00000001/qq.dac.o.00000001",
        "remap/qq/00000/00000000/qq.dac.00000003.00000000",
        "remap/qq/00000/00000001/qq.dac.00000003.00000001",
        "slice/params.dac",
        "slice/qqz.dac.00000004.00000001",
        "slice/slice.dac",
        "tau/qq.dac.00000004",
        "time/tau/t.dac.00000004",
    ]
