"""pyR2D2.Sync の zarr の部分取り出し (remap_zarr, tau_zarr, slice_zarr) の試験。

ssh を「接続先を捨てて手元のシェルで動かす」スクリプトに差し替え、リモートのスクリプトを
実際に手元の Python で動かす。リモートの zarr は .zarr.zip とディレクトリの両方を試す。
"""

import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

zarr = pytest.importorskip("zarr")

from pyR2D2.data_io import zarr_util
from pyR2D2.sync.sync import Sync

CASEID = "d001"
SHAPE = (10, 6, 8)
CHUNKS = (4, 3, 8)
KEYS = ["ro", "vx", "se"]


@pytest.fixture
def fake_ssh(tmp_path):
    """``fake_ssh host command`` で、host を捨てて command を手元の sh で動かす。"""
    path = tmp_path / "fake_ssh"
    path.write_text('#!/bin/sh\nshift\nexec sh -c "$*"\n')
    path.chmod(0o755)
    return str(path)


def field(key, n):
    base = {"ro": 1.0, "vx": 2.0, "se": 3.0}[key]
    return (base + n + np.arange(np.prod(SHAPE)).reshape(SHAPE) / 1000).astype(np.float32)


def make_remote(root, steps, zipped):
    data = root / CASEID / "data"
    for n in steps:
        path = data / "remap" / "qq" / "zarr" / f"qq.{n:08d}.zarr"
        vars_dict = {k: field(k, n) for k in KEYS}
        vars_dict.update(x=np.linspace(0, 1, SHAPE[0]), y=np.arange(SHAPE[1]), z=np.arange(SHAPE[2]))
        zarr_util.save(path, vars_dict, {"n": n}, chunks3d=CHUNKS)
        if zipped:
            zarr_util.zip_zarr(path, remove_original=True)
    tau = data / "tau" / "zarr" / f"qq.{steps[0]:08d}.zarr"
    zarr_util.save(tau, {k: field(k, 0)[0] for k in KEYS}, {})
    sl = data / "slice" / "zarr" / "x" / f"qq.{steps[0]:08d}.zarr"
    zarr_util.save(sl, {k + "_yin": field(k, 0)[:2] for k in KEYS}, {})
    return data


def make_sync(tmp_path):
    datadir = tmp_path / "local" / CASEID / "data"
    datadir.mkdir(parents=True)
    x = np.linspace(0, 1, SHAPE[0])
    data = SimpleNamespace(datadir=datadir, p=SimpleNamespace(), x=x, y=np.arange(SHAPE[1]), z=np.arange(SHAPE[2]))
    return Sync(data), datadir


@pytest.mark.parametrize("zipped", [True, False])
def test_remap_zarr_series_region_and_keys(tmp_path, fake_ssh, zipped):
    remote_root = tmp_path / "remote" / "run"
    make_remote(remote_root, [1, 2, 3], zipped)
    sync, datadir = make_sync(tmp_path)

    written = sync.remap_zarr(
        "host", [1, 2, 3, 4], keys=["ro", "vx"], i0=5, i1=6, j0=1, j1=4,
        remote_root=remote_root, remote_python=sys.executable, ssh=fake_ssh,
    )
    assert [p.name for p in written] == [f"qq.{n:08d}.zarr" for n in [1, 2, 3]]

    for n in [1, 2, 3]:
        path = datadir / "remap" / "qq" / "zarr" / f"qq.{n:08d}.zarr"
        qq, params = zarr_util.load(path, names="all", with_attrs=True)
        assert params == {"n": n}
        assert set(qq) == {"ro", "vx", "x", "y", "z"}
        np.testing.assert_array_equal(qq["x"][:], np.linspace(0, 1, SHAPE[0]).astype(np.float32))
        for key in ["ro", "vx"]:
            assert qq[key].shape == SHAPE
            np.testing.assert_array_equal(qq[key][5, 1:4], field(key, n)[5, 1:4])
            assert np.isnan(qq[key][4]).all() and np.isnan(qq[key][5, 0]).all()
        # 添字を変えずに Data の qx/qr と同じ読み方ができる
        layer = zarr_util.load(path, names=["ro"], i0=5, i1=6, j0=1, j1=4)["ro"]
        np.testing.assert_array_equal(layer, field("ro", n)[5:6, 1:4])


def test_remap_zarr_skips_done_and_merges(tmp_path, fake_ssh, capsys):
    remote_root = tmp_path / "remote" / "run"
    make_remote(remote_root, [1], zipped=True)
    sync, datadir = make_sync(tmp_path)
    kw = dict(remote_root=remote_root, remote_python=sys.executable, ssh=fake_ssh)

    assert len(sync.remap_zarr("host", 1, keys=["ro"], x0=0.5, x1=0.6, **kw)) == 1
    # 同じ範囲 (の一部) はもう来ない
    assert sync.remap_zarr("host", 1, keys=["ro"], i0=4, i1=5, **kw) == []
    assert "already here" in capsys.readouterr().out
    # 別の層と変数は同じ手元の zarr に足される
    assert len(sync.remap_zarr("host", 1, keys=["se"], i0=0, i1=1, **kw)) == 1

    path = datadir / "remap" / "qq" / "zarr" / "qq.00000001.zarr"
    root = zarr_util.open_zarr_group(path)
    i = np.argmin(abs(np.linspace(0, 1, SHAPE[0]) - 0.5))
    np.testing.assert_array_equal(root["ro"][i : i + 2], field("ro", 1)[i : i + 2])
    np.testing.assert_array_equal(root["se"][0], field("se", 1)[0])
    assert np.isnan(root["se"][i]).all()
    assert len(root.attrs["pyR2D2_partial"]) == 2


def test_tau_and_slice_zarr(tmp_path, fake_ssh):
    remote_root = tmp_path / "remote" / "run"
    make_remote(remote_root, [7], zipped=False)
    sync, datadir = make_sync(tmp_path)
    kw = dict(remote_root=remote_root, remote_python=sys.executable, ssh=fake_ssh)

    sync.tau_zarr("host", 7, keys=["vx"], j0=2, j1=3, **kw)
    tau = zarr_util.open_zarr_group(datadir / "tau" / "zarr" / "qq.00000007.zarr")
    assert list(tau.array_keys()) == ["vx"]
    np.testing.assert_array_equal(tau["vx"][2], field("vx", 0)[0, 2])

    sync.slice_zarr("host", 7, "x", keys=["se"], n_slice=1, **kw)
    sl = zarr_util.open_zarr_group(datadir / "slice" / "zarr" / "x" / "qq.00000007.zarr")
    assert list(sl.array_keys()) == ["se_yin"]
    np.testing.assert_array_equal(sl["se_yin"][1], field("se", 0)[1])
    assert np.isnan(sl["se_yin"][0]).all()


def test_remote_failure_is_reported(tmp_path, fake_ssh):
    remote_root = tmp_path / "remote" / "run"
    make_remote(remote_root, [1], zipped=True)
    sync, datadir = make_sync(tmp_path)
    with pytest.raises(RuntimeError):
        sync.remap_zarr(
            "host", 1, keys=["nosuch"], remote_root=remote_root,
            remote_python=sys.executable, ssh=fake_ssh,
        )
    assert not (datadir / "remap" / "qq" / "zarr" / "qq.00000001.zarr").exists()
