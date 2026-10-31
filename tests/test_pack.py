"""r2d2plus-pack (pyR2D2.pack) の試験。

R2D2plus の見本 (tests/data/r2d2plus_fixtures、compressed と legacy) を 2 つ展開し、片方をまとめて元を消す。
wall_limit で止まったランとして扱うと最後のステップは元が残るので、まとめた方は「zip のステップと
元のファイルのステップが混ざったラン」になる。
これを読み手 (qf, qx, qz, qm, qr, qs, qt) で読み、まとめていない方とビット単位で同じことを確かめる。
"""

import tarfile
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("zarr")

import pyR2D2
from pyR2D2.pack import Packer, parse_steps
from pyR2D2.pack.archive import Archiver
from pyR2D2.pack.cli import archive_main
from pyR2D2.pack.cli import main as pack_main

try:
    import tomllib
except ImportError:  # pragma: no cover
    tomllib = None

FIXTURE = sorted((Path(__file__).parent / "data" / "r2d2plus_fixtures").glob("r2d2plus-*.tar.xz"))[-1]
TAUS = ["", "01", "001"]


def extract(dest):
    with tarfile.open(FIXTURE) as tar:
        tar.extractall(dest, filter="data")
    return dest


@pytest.fixture(params=["compressed", "legacy"])
def runs(request, tmp_path):
    if request.param == "compressed":
        # 新出力形式の zstd は numcodecs か zstandard で解く
        try:
            import numcodecs  # noqa: F401
        except ImportError:
            pytest.importorskip("zstandard")
    orig = extract(tmp_path / "orig") / request.param / "data"
    packed = extract(tmp_path / "packed") / request.param / "data"
    return request.param, orig, packed


def report(datadir):
    with open(datadir / "param" / "pack_report.toml", "rb") as f:
        return tomllib.load(f)["pack"]


def stored_qq(d):
    if d.p.output_format == "compressed":
        return list(d.p.format_info["variables"])
    return d.remap_keys + d.remap_keys_add


def same(a, b, what):
    a, b = np.asarray(a), np.asarray(b)
    assert a.shape == b.shape, what
    assert a.dtype.itemsize == 4 and b.dtype.itemsize == 4, what
    np.testing.assert_array_equal(
        a.astype("=f4").view(np.uint32), b.astype("=f4").view(np.uint32), err_msg=what
    )


def test_parse_steps():
    assert parse_steps(None) is None
    assert parse_steps("1,3-5") == {1, 3, 4, 5}


def set_status(datadir, status):
    summary = datadir / "param" / "run_summary.toml"
    text = summary.read_text()
    old = [line for line in text.splitlines() if line.startswith("status = ")][0]
    summary.write_text(text.replace(old, f'status = "{status}"', 1))


def files(root):
    return sorted(str(p.relative_to(root)) for p in Path(root).rglob("*"))


def test_pack_delete_and_read_mixed(runs):
    fmt, orig, packed = runs
    set_status(packed, "wall_limit")
    assert Packer(packed).run(delete_original=True)

    entries = report(packed)
    assert {e["result"] for e in entries} == {"packed"}
    last = {(e["kind"], e.get("direc")): max(f["n"] for f in entries if (f["kind"], f.get("direc")) == (e["kind"], e.get("direc"))) for e in entries}
    for e in entries:
        is_last = e["n"] == last[(e["kind"], e.get("direc"))]
        assert e["deleted_files"] == (0 if is_last else e["raw_files"]), e
        assert (packed / e["zip"]).is_file()
    # 3D は元に書かれた変数だけ (新出力形式では pr, te, op を入れない)
    qf_keys = [e for e in entries if e["kind"] == "qf"][0]["keys"]
    if fmt == "compressed":
        assert qf_keys == ["ro", "vx", "vy", "vz", "bx", "by", "bz", "se"]
    # 作業用の .packing は残らない
    assert not list(packed.rglob("*.packing*"))
    # 消したのは照合したステップの元のファイルだけ
    assert not list((packed / "remap" / "qq").glob("*/*/qq.*.00000000.*"))
    assert list((packed / "remap" / "qq").glob("*/*/qq.*.00000002.*"))

    do, dp = pyR2D2.Data(str(orig)), pyR2D2.Data(str(packed))
    keys = stored_qq(do)
    for n in range(3):
        do.qf.read(n)
        dp.qf.read(n)
        for k in keys:
            same(getattr(dp.qf, k), getattr(do.qf, k), f"qf {n} {k}")
        for reader, args in [
            ("qx", (do.x[3], n)),
            ("qz", (do.z[1], n)),
            ("qm", (0, n)),
        ]:
            getattr(do, reader).read(*args)
            getattr(dp, reader).read(*args)
            for k in keys:
                same(getattr(getattr(dp, reader), k), getattr(getattr(do, reader), k), f"{reader} {n} {k}")
        kw = dict(x0=do.x[5], x1=do.x[20], y0=do.y[1], y1=do.y[2])
        do.qr.read(n, keys, **kw)
        dp.qr.read(n, keys, **kw)
        for k in keys:
            same(getattr(dp.qr, k), getattr(do.qr, k), f"qr {n} {k}")

        do.qt.read(n)
        dp.qt.read(n)
        tau_keys = [e for e in entries if e["kind"] == "qt"][0]["keys"]
        for k in tau_keys:
            same(getattr(dp.qt, k), getattr(do.qt, k), f"qt {n} {k}")
        for t in TAUS:
            np.testing.assert_array_equal(getattr(dp.qt, "he" + t), getattr(do.qt, "he" + t))

        slice_keys = [e for e in entries if e["kind"] == "qs"][0]["keys"]
        for direc in "xyz":
            for i in range(len(do.qs.xyz_slice_select(direc))):
                do.qs.read(i, direc, n)
                dp.qs.read(i, direc, n)
                for k in slice_keys:
                    same(getattr(dp.qs, k), getattr(do.qs, k), f"qs {direc} {i} {n} {k}")


def test_normal_end_deletes_last_step_too(runs):
    fmt, orig, packed = runs
    assert Packer(packed).run(kinds=["qf", "qt"], delete_original=True)
    assert all(e["deleted_files"] == e["raw_files"] for e in report(packed))
    assert not list((packed / "remap" / "qq").glob("*/*/qq.*"))
    d = pyR2D2.Data(str(packed))
    d.qf.read(2)
    assert d.qf.ro.shape == (d.ix, d.jx, d.kx)


def test_out_leaves_run_untouched(runs, tmp_path):
    fmt, orig, packed = runs
    before = files(packed)
    out = tmp_path / "out"
    assert pack_main([str(packed), "--kinds", "qf,qs", "--steps", "1", "--out", str(out)]) == 0
    assert files(packed) == before
    assert (out / "data" / "remap" / "qq" / "zarr" / "qq.00000001.zarr.zip").is_file()
    assert (out / "data" / "slice" / "zarr" / "x" / "qq.00000001.zarr.zip").is_file()
    assert [e["zip"] for e in report(out / "data")][0] == "remap/qq/zarr/qq.00000001.zarr.zip"
    with pytest.raises(ValueError):
        Packer(packed, out=out).run(delete_original=True)
    with pytest.raises(ValueError):
        Packer(packed, out=packed.parent / "x")


def test_archive(runs, tmp_path):
    fmt, orig, packed = runs
    import tarfile, zipfile
    (packed / "param" / "source").mkdir(exist_ok=True)
    (packed / "param" / "source" / "src.tar.zst").write_text("x")
    out = tmp_path / "archive"
    # まとめていない出力が残っていれば止める
    assert archive_main([str(packed.parent), str(out)]) == 1
    assert not out.exists()
    assert Packer(packed).run(delete_original=True)
    before = files(packed.parent)
    assert archive_main([str(packed.parent), str(out)]) == 0
    assert files(packed.parent) == before  # 元のランには触らない

    slots = sorted(p.name for p in (packed / "restart").iterdir())
    assert sorted(p.name for p in (out / "data" / "restart").iterdir()) == [s + ".tar" for s in slots]
    with tarfile.open(out / "data" / "restart" / f"{slots[0]}.tar") as tf:
        assert any(m.name.endswith("meta.toml") for m in tf)
    with zipfile.ZipFile(out / "data" / "time.zip") as zf:
        assert zf.testzip() is None
        assert sorted(zf.namelist()) == sorted(
            str(p.relative_to(packed / "time")) for p in (packed / "time").rglob("*") if p.is_file()
        )
    with tarfile.open(out / "rest.tar") as tf:
        names = {m.name for m in tf}
    assert f"{fmt}/data/param/source/src.tar.zst" not in names  # run の相対で入る
    assert "data/param/source/src.tar.zst" in names and "data/param/params.dac" in names
    assert not any(n.startswith(("data/restart", "data/time", "data/remap/qq", "data/tau/", "data/slice/zarr")) for n in names)
    assert not any(n.startswith("data/slice/qq") for n in names)
    assert "data/slice/params.dac" in names
    with open(out / "archive_report.toml", "rb") as f:
        rep = tomllib.load(f)
    assert "remap/qq/zarr/qq.00000000.zarr.zip" in rep["zips"]
    assert {p["part"] for p in rep["part"]} == {"restart", "time", "rest"}

    set_status(packed, "running")
    with pytest.raises(RuntimeError):
        Archiver(packed.parent, tmp_path / "again").run_all()
    with pytest.raises(ValueError):
        Archiver(packed.parent, packed / "x")


def test_verify_existing_then_delete(runs, capsys):
    fmt, orig, packed = runs
    # まず消さずにまとめる。次に消すと、前からある zip を元と照合してから消す
    assert pack_main([str(packed.parent), "--kinds", "qf", "--workers", "1"]) == 0
    assert all(e["kept_because"] == "not requested" for e in report(packed))
    assert pack_main([str(packed), "--kinds", "qf", "--delete-original", "--dry-run"]) == 0
    assert "verify" in capsys.readouterr().out
    assert len(report(packed)) == 3  # dry-run は記録しない
    assert pack_main([str(packed), "--kinds", "qf", "--delete-original", "--steps", "0-1"]) == 0
    second = report(packed)[3:]
    assert [e["result"] for e in second] == ["verified", "verified"]
    assert [e["deleted_files"] for e in second] == [1, 1]


def test_running_run_packs_finished_steps(runs):
    """走っているランは、後の時刻のファイルがあるステップ (書き終わったもの) だけをまとめて元を消す。"""
    fmt, orig, packed = runs
    set_status(packed, "running")
    assert Packer(packed).run(kinds=["qf", "qt"], delete_original=True)
    done = {(e["kind"], e["n"]) for e in report(packed)}
    assert done == {("qf", 0), ("qf", 1), ("qt", 0), ("qt", 1)}
    assert all(e["deleted_files"] == e["raw_files"] for e in report(packed))
    # いちばん新しいステップ (後の時刻が無い) はまとめず、元を残す
    assert not (packed / "remap" / "qq" / "zarr" / "qq.00000002.zarr.zip").exists()
    assert list((packed / "remap" / "qq").glob("*/*/qq.*.00000002.*"))
    # restart と time には触らない
    assert list((packed / "time" / "mhd").glob("t.dac.*")) and list((packed / "restart").iterdir())

    # 時刻が time.zip にしかなくても見分ける。次の出力が出れば n=2 も書き終わった扱いになる
    import shutil, zipfile
    with zipfile.ZipFile(packed / "time.zip", "w") as zf:
        for f in sorted((packed / "time").rglob("*")):
            if f.is_file():
                zf.write(f, f.relative_to(packed))
        zf.writestr("time/mhd/t.dac.00000003", b"\0" * 8)
    shutil.rmtree(packed / "time")
    assert Packer(packed).run(kinds=["qf"], delete_original=True)
    assert (packed / "remap" / "qq" / "zarr" / "qq.00000002.zarr.zip").exists()


def test_mismatch_keeps_originals(runs):
    fmt, orig, packed = runs
    assert Packer(packed).run(kinds=["qf"])
    # zip を別のステップのものにすり替えると、照合に失敗して元は消さない
    z = packed / "remap" / "qq" / "zarr"
    (z / "qq.00000000.zarr.zip").write_bytes((z / "qq.00000001.zarr.zip").read_bytes())
    assert not Packer(packed).run(kinds=["qf"], steps={0}, delete_original=True)
    e = report(packed)[-1]
    assert e["result"] == "failed" and e["deleted_files"] == 0
    assert list((packed / "remap" / "qq").glob("*/*/qq.*.00000000.*"))


def test_incomplete_step_is_skipped(runs):
    fmt, orig, packed = runs
    for f in (packed / "slice").glob("qqx*.00000001.00000001"):
        f.unlink()
    Packer(packed).run(kinds=["qs"], delete_original=True)
    assert not (packed / "slice" / "zarr" / "x" / "qq.00000001.zarr.zip").exists()
    assert list((packed / "slice").glob("qqx*.00000001.00000002"))


def test_restore(runs, tmp_path):
    """archive の restart の tar を戻し、archive_report.toml の一覧と照合して展開する。"""
    import filecmp, shutil
    from pyR2D2.pack.cli import restore_main

    fmt, orig, packed = runs
    assert Packer(packed).run(delete_original=True)
    out = tmp_path / "archive"
    assert archive_main([str(packed.parent), str(out), "--parts", "restart"]) == 0
    slot = sorted(p.name for p in (packed / "restart").iterdir())[0]
    original = tmp_path / "slot_copy"
    shutil.copytree(packed / "restart" / slot, original, symlinks=True)

    # すでにあれば戻さない
    assert restore_main([str(packed.parent), "--slot", slot, "--from", str(out)]) == 1
    shutil.rmtree(packed / "restart" / slot)
    assert restore_main([str(packed.parent), "--slot", slot, "--from", str(out), "--dry-run"]) == 0
    assert not (packed / "restart" / slot).exists()
    assert restore_main([str(packed.parent), "--slot", slot, "--from", str(out)]) == 0
    cmp = filecmp.dircmp(original, packed / "restart" / slot)
    assert not cmp.diff_files and not cmp.left_only and not cmp.right_only
    assert not list((packed / "restart").glob(".restore-*"))
    with open(packed / "param" / "restore_report.toml", "rb") as f:
        assert tomllib.load(f)["restore"][0]["slot"] == slot
    assert restore_main([str(packed.parent), "--slot", "nosuch", "--from", str(out)]) == 1

    # 壊れた tar (SHA-256 が違う) は展開しない
    shutil.rmtree(packed / "restart" / slot)
    tar = out / "data" / "restart" / f"{slot}.tar"
    data = bytearray(tar.read_bytes())
    data[-1] ^= 1
    tar.write_bytes(bytes(data))
    assert restore_main([str(packed.parent), "--slot", slot, "--from", str(out)]) == 1
    assert not (packed / "restart" / slot).exists()


def test_restore_from_gfarm_uses_gfpcopy(runs, tmp_path):
    """gfarm: の場所からは gfpcopy で持ってくる (ここでは手元のパスに読み替える偽物)。"""
    import shutil
    from pyR2D2.pack.restore import Restorer

    fmt, orig, packed = runs
    assert Packer(packed).run(delete_original=True)
    out = tmp_path / "archive"
    assert archive_main([str(packed.parent), str(out), "--parts", "restart"]) == 0
    slot = sorted(p.name for p in (packed / "restart").iterdir())[-1]
    shutil.rmtree(packed / "restart" / slot)
    fake = tmp_path / "gfpcopy"
    fake.write_text(f'#!/bin/sh\nshift\nsrc=$(echo "$1" | sed "s|^gfarm:///place|{out}|")\ncp "$src" "$2"/\n')
    fake.chmod(0o755)
    r = Restorer(packed.parent, source="gfarm:///place", gfpcopy=str(fake), log=lambda *m: None)
    assert r.run_slot(slot) == packed / "restart" / slot


def test_restore_default_source():
    from pyR2D2.pack.restore import default_source

    assert default_source("/vol0600/data/hp200137/u10094/R2D2plus/run/d00012") == (
        "gfarm:///home/hp200137/hpci000196/fugaku/R2D2plus/run/d00012"
    )
    assert default_source("/x/run/d001", project="lowmass") == "gfarm:///home/hp200137/hpci000196/fugaku/lowmass/run/d001"
    with pytest.raises(ValueError):
        default_source("/x/somewhere/d001")
