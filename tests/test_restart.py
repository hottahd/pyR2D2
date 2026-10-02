"""R2D2plus の restart の読み手 (pyR2D2.read_restart) の試験。

合成した restart で、v1 / v2・ランク分割の組み立て・YinYang・圧縮・照合・最新の slot の選び方を確かめる。
実物の restart での照合 (1 ランクと鉛直 2 分割が完全一致すること) は、R2D2plus の verify_run があるときだけ走る。
"""

import os
from pathlib import Path

import numpy as np
import pytest

import pyR2D2
from pyR2D2.data_io import compressed
from pyR2D2.data_io.restart import RESTART_VARIABLES


def global_state(shape, seed=0):
    rng = np.random.default_rng(seed)
    return {n: rng.standard_normal(shape) * (i + 1) for i, n in enumerate(RESTART_VARIABLES)}


def write_restart(slot, state, mpi=(1, 1, 1), version=2, step=10, time=1.5, compress=False,
                  yinyang=False, extras=None, checksum=True):
    """R2D2plus と同じ並びで restart を書く (試験用の書き手)。"""
    slot = Path(slot)
    slot.mkdir(parents=True)
    panels = {"yin": state[0], "yang": state[1]} if yinyang else {None: state}
    first = next(iter(panels.values()))
    nx, ny, nz = first["ro"].shape
    sx, sy, sz = mpi
    lx, ly, lz = nx // sx, ny // sy, nz // sz
    gen = ".generation-TEST01"
    payload = slot / gen if version == 2 else slot
    digests = []
    for panel, st in panels.items():
        d = payload / panel if panel else payload
        d.mkdir(parents=True, exist_ok=True)
        for rank in range(sx * sy * sz):
            ix, rem = divmod(rank, sy * sz)
            iy, iz = divmod(rem, sz)
            local = np.stack([st[n][ix * lx:(ix + 1) * lx, iy * ly:(iy + 1) * ly, iz * lz:(iz + 1) * lz]
                              for n in RESTART_VARIABLES])
            f = d / f"state.{rank:08d}.bin"
            if compress:
                compressed.write_r2d2plus_z(f, "restart", [("state", local.reshape(-1), "C")])
            else:
                f.write_bytes(local.astype("<f8").tobytes())
            digests.append(compressed._xxh3_hex(f.read_bytes()))
    lines = ["[restart]", f"format_version = {version}"]
    if version == 2:
        lines.append(f'generation = "{gen}"')
    lines += [f"step = {step}", f"time = {time!r}", "", "[grid]", f"nx = {nx}", f"ny = {ny}", f"nz = {nz}",
              "", "[mpi]", f"size_x = {sx}", f"size_y = {sy}", f"size_z = {sz}"]
    if yinyang:
        lines.append("yinyang = true")
    if checksum and all(digests):
        lines += ["", "[checksum]", 'algorithm = "xxh3_64"', "values = [" + ", ".join(f'"{d}"' for d in digests) + "]"]
    for name, data in (extras or {}).items():
        (payload / f"extra.{name}.bin").write_bytes(data)
        lines += ["", f"[extras.{name}]", "shared = true", 'values = ["0"]']
    (slot / "meta.toml").write_text("\n".join(lines) + "\n")


@pytest.mark.parametrize("mpi", [(1, 1, 1), (2, 1, 1), (1, 2, 1), (1, 1, 2), (2, 3, 2)])
@pytest.mark.parametrize("version", [1, 2])
def test_assembles_ranks(tmp_path, mpi, version):
    st = global_state((4, 6, 8))
    write_restart(tmp_path / "data" / "restart" / "e", st, mpi=mpi, version=version)
    r = pyR2D2.read_restart(tmp_path)
    assert r.shape == (4, 6, 8) and r.step == 10 and r.time == 1.5
    for n in RESTART_VARIABLES:
        np.testing.assert_array_equal(r.state[n], st[n])


def test_compressed_state(tmp_path):
    pytest.importorskip("xxhash")
    pytest.importorskip("zstandard")
    st = global_state((4, 4, 6), seed=1)
    write_restart(tmp_path / "e", st, mpi=(2, 1, 1), compress=True)
    r = pyR2D2.read_restart(tmp_path / "e", verify=True)
    for n in RESTART_VARIABLES:
        np.testing.assert_array_equal(r.state[n], st[n])


def test_yinyang_and_extras(tmp_path):
    yin, yang = global_state((2, 6, 4), seed=2), global_state((2, 6, 4), seed=3)
    write_restart(tmp_path / "o", (yin, yang), mpi=(1, 2, 1), yinyang=True, extras={"afc": b"\x01\x02\x03"})
    r = pyR2D2.read_restart(tmp_path)
    assert r.yinyang and set(r.state) == {"yin", "yang"}
    np.testing.assert_array_equal(r.state["yin"]["se"], yin["se"])
    np.testing.assert_array_equal(r.state["yang"]["bz"], yang["bz"])
    assert r.extra_bytes("afc") == b"\x01\x02\x03"


def test_latest_slot_and_listing(tmp_path):
    root = tmp_path / "data" / "restart"
    write_restart(root / "e", global_state((2, 2, 2)), step=30, time=3.0)
    write_restart(root / "o", global_state((2, 2, 2)), step=20, time=2.0)
    write_restart(root / "00000001", global_state((2, 2, 2)), step=10, time=1.0)
    assert pyR2D2.restart_slots(tmp_path) == {"00000001": (10, 1.0), "e": (30, 3.0), "o": (20, 2.0)}
    assert pyR2D2.read_restart(tmp_path).path.name == "e"
    assert pyR2D2.read_restart(tmp_path, slot="o").step == 20
    with pytest.raises(FileNotFoundError):
        pyR2D2.read_restart(tmp_path, slot="00000099")


def test_checksum_mismatch_is_detected(tmp_path):
    pytest.importorskip("xxhash")
    write_restart(tmp_path / "e", global_state((2, 2, 4)))
    f = next((tmp_path / "e").rglob("state.00000000.bin"))
    raw = bytearray(f.read_bytes())
    raw[0] ^= 0xFF
    f.write_bytes(bytes(raw))
    with pytest.raises(ValueError, match="checksum"):
        pyR2D2.read_restart(tmp_path / "e")
    pyR2D2.read_restart(tmp_path / "e", verify=False)  # 照合を切れば読める


def test_wrong_size_and_version_are_rejected(tmp_path):
    write_restart(tmp_path / "e", global_state((2, 2, 4)), checksum=False)
    f = next((tmp_path / "e").rglob("state.00000000.bin"))
    f.write_bytes(f.read_bytes()[:-8])
    with pytest.raises(ValueError, match="expected"):
        pyR2D2.read_restart(tmp_path / "e")
    meta = tmp_path / "e" / "meta.toml"
    meta.write_text(meta.read_text().replace("format_version = 2", "format_version = 3"))
    with pytest.raises(ValueError, match="format_version"):
        pyR2D2.read_restart(tmp_path / "e")


REAL = Path(os.environ.get("R2D2PLUS_VERIFY_RUN", "/scr/a000/c0234hotta/Repository/R2D2plus/verify_run"))


@pytest.mark.skipif(not (REAL / "artdif_upper_seam_20260929" / "fix40").is_dir(), reason="R2D2plus の verify_run が無い")
def test_real_decomposition_independence():
    """R2D2plus DEC-586 の修正後、1 ランクと鉛直 2 分割の 40 step 目は完全一致する (読み手の組み立ての検証)。"""
    base = REAL / "artdif_upper_seam_20260929" / "fix40"
    a = pyR2D2.read_restart(base / "r1", slot="e")
    b = pyR2D2.read_restart(base / "r2", slot="e")
    assert a.meta["mpi"]["size_z"] == 1 and b.meta["mpi"]["size_z"] == 2
    for n in RESTART_VARIABLES:
        np.testing.assert_array_equal(a.state[n], b.state[n])
