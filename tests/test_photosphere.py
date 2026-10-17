"""R2D2plus の光球スラブ出力の読み手 (pyR2D2.Photosphere / read_photosphere) の試験。

R2D2plus の PhotosphereWriter と同じ並びで書いた合成の出力で、版 1 (6 量固定、``units`` の無い旧 format.toml) と
版 2 (磁場を含む ``fields``) の読み分け、ランク分割の組み立て (鉛直の範囲に掛からないランクを含む)、
量と範囲の部分読み、壊れた出力の拒否を確かめる。
"""

from pathlib import Path

import numpy as np
import pytest

import pyR2D2
from pyR2D2.data_io.photosphere import PHOTOSPHERE_FIELDS

# 全体の格子 (鉛直 nz_total、水平 nx, ny) とランク分割。スラブは鉛直の [K0, K1) だけ。
NZ_TOTAL, NX, NY = 8, 6, 10
MPI = (2, 2, 2)  # (x, y, z)
K0, K1 = 1, 4    # 下の z ランク (0..3) と上の z ランク (4..7) のうち、上のランクは掛からない
MPI_Z_LOCAL = NZ_TOTAL // MPI[2]


def global_values(names, shape):
    """量ごとに違う、符号の混じった半端な値 (float32 に丸めた後の値)。"""
    nz, nx, ny = shape
    k, i, j = np.meshgrid(np.arange(nz), np.arange(nx), np.arange(ny), indexing="ij")
    out = {}
    for n, name in enumerate(names):
        p = PHOTOSPHERE_FIELDS.index(name)
        out[name] = (0.37 * (p + 1) * 1.0e3 - 41.3 * k + 7.9 * (p + 2) * i - 2.3 * j * (p - 3.5)).astype(np.float32)
    return out


def write_output(directory, names, version, header_comment, units=True, snapshots=((0, 10, 12.5),)):
    """R2D2plus の PhotosphereWriter と同じ並びで書く (試験用の書き手)。"""
    directory = Path(directory)
    directory.mkdir(parents=True)
    shape = (K1 - K0, NX, NY)
    np.linspace(0.0, 1.2e8, NX).astype("<f8").tofile(directory / "x.f64")
    np.linspace(0.0, 2.0e8, NY).astype("<f8").tofile(directory / "y.f64")
    np.linspace(6.95e10, 6.96e10, NZ_TOTAL)[K0:K1].astype("<f8").tofile(directory / "z.f64")
    sx, sy, sz = MPI
    lx, ly = NX // sx, NY // sy
    fields = ", ".join(f'"{n}"' for n in names)
    text = (f"# {header_comment}\n\n[format]\nversion = {version}\ndtype = \"float32\"\nbyte_order = \"little\"\n"
            "layout = \"field-major, C-order within each field\"\naxes = [\"z\", \"x\", \"y\"]\n"
            f"fields = [{fields}]\n")
    if units:
        text += "units = [" + ", ".join(f'"{pyR2D2.data_io.photosphere.PHOTOSPHERE_UNITS[n]}"' for n in names) + "]\n"
    text += (f"coordinate_dtype = \"float64\"\ncomplete_marker = \"meta.toml\"\n\n[grid]\n"
             f"shape = [{shape[0]}, {NX}, {NY}]\nglobal_k_begin = {K0}\nglobal_k_end = {K1}\n"
             f"x = \"x.f64\"\ny = \"y.f64\"\nz = \"z.f64\"\n\n[mpi]\nsize = [{sx}, {sy}, {sz}]\n"
             "rank_files = \"state.{rank:04d}.f32\"\n")
    ranks = []
    for rank in range(sx * sy * sz):
        cx, rem = divmod(rank, sy * sz)
        cy, cz = divmod(rem, sz)
        lo = max(K0, cz * MPI_Z_LOCAL)
        hi = min(K1, (cz + 1) * MPI_Z_LOCAL)
        nz = max(0, hi - lo)
        ranks.append((rank, (lo - K0, cx * lx, cy * ly), (nz, lx, ly)))
        text += (f"\n[[rank]]\nid = {rank}\ncoords = [{cx}, {cy}, {cz}]\n"
                 f"start = [{lo - K0}, {cx * lx}, {cy * ly}]\nshape = [{nz}, {lx}, {ly}]\n")
    (directory / "format.toml").write_text(text)

    expected = {}
    for index, step, time in snapshots:
        values = global_values(names, shape)
        values = {n: v + np.float32(0.11 * index) for n, v in values.items()}
        expected[index] = values
        snap = directory / f"{index:08d}"
        snap.mkdir()
        for rank, start, rshape in ranks:
            sl = tuple(slice(s, s + n) for s, n in zip(start, rshape))
            with open(snap / f"state.{rank:04d}.f32", "wb") as f:
                for n in names:
                    f.write(np.ascontiguousarray(values[n][sl]).astype("<f4").tobytes())
        (snap / "meta.toml").write_text(
            "# Presence of this file marks a complete snapshot.\n\n[snapshot]\n"
            f"format_version = {version}\nindex = {index}\nstep = {step}\ntime = {time!r}\n")
    return expected


def test_reads_version1_six_fields_written_by_the_old_writer(tmp_path):
    # 2026-10-17 までの R2D2plus が書いた format.toml そのまま (units が無い、見出しの DEC 番号も旧い)
    root = tmp_path / "data" / "photosphere"
    expected = write_output(root, PHOTOSPHERE_FIELDS[:6], 1, "R2D2+ photosphere snapshot format (DEC-443)",
                            units=False, snapshots=((0, 10, 12.5), (3, 47, 27.25)))
    ph = pyR2D2.Photosphere(tmp_path)  # ランのディレクトリからでも見つける
    assert ph.version == 1
    assert ph.fields == PHOTOSPHERE_FIELDS[:6]
    assert ph.units["pressure"] == "dyn/cm^2"
    assert ph.shape == (K1 - K0, NX, NY)
    assert ph.snapshots() == [0, 3]

    s = pyR2D2.read_photosphere(root)  # 既定は最後の番号
    assert (s.index, s.step, s.time) == (3, 47, 27.25)
    for name in PHOTOSPHERE_FIELDS[:6]:
        np.testing.assert_array_equal(s[name], expected[3][name])
    np.testing.assert_array_equal(s.z, np.linspace(6.95e10, 6.96e10, NZ_TOTAL)[K0:K1])
    assert s.x.shape == (NX,) and s.y.shape == (NY,)


def test_reads_version2_with_magnetic_field_and_subsets(tmp_path):
    names = PHOTOSPHERE_FIELDS  # 版 1 の 6 量 + 磁場 3 成分
    expected = write_output(tmp_path / "ph", names, 2,
                            "R2D2+ photosphere snapshot format (docs/formats/photosphere.md, DEC-641)")
    ph = pyR2D2.Photosphere(tmp_path / "ph")
    assert ph.version == 2
    assert ph.units["magnetic_vertical"] == "G"
    full = ph.read(0)
    for name in names:
        np.testing.assert_array_equal(full[name], expected[0][name])

    # 量と範囲の部分読み (ランクの境界をまたぐ範囲)
    part = ph.read(0, fields=["magnetic_horizontal_y", "velocity_vertical"], z=slice(-2, None), x=(2, 5), y=slice(3, 9))
    assert set(part.fields) == {"magnetic_horizontal_y", "velocity_vertical"}
    assert part.units == {"magnetic_horizontal_y": "G", "velocity_vertical": "cm/s"}
    for name in part.fields:
        np.testing.assert_array_equal(part[name], expected[0][name][-2:, 2:5, 3:9])
    np.testing.assert_array_equal(part.x, ph.x[2:5])
    np.testing.assert_array_equal(part.z, ph.z[-2:])


def test_reads_a_subset_written_without_thermodynamics(tmp_path):
    names = ("velocity_vertical", "magnetic_vertical", "magnetic_horizontal_x")
    expected = write_output(tmp_path / "ph", names, 2, "subset")
    s = pyR2D2.read_photosphere(tmp_path / "ph", 0)
    assert tuple(s.fields) == names
    for name in names:
        np.testing.assert_array_equal(s[name], expected[0][name])
    with pytest.raises(KeyError):
        pyR2D2.read_photosphere(tmp_path / "ph", 0, fields=["temperature"])


def test_rejects_unknown_version_truncated_files_and_incomplete_snapshots(tmp_path):
    root = tmp_path / "ph"
    write_output(root, PHOTOSPHERE_FIELDS[:6], 1, "v1", units=False)
    ph = pyR2D2.Photosphere(root)
    # 書きかけ (meta.toml が無い) は番号の一覧に出ず、読むと止まる
    (root / "00000001").mkdir()
    assert ph.snapshots() == [0]
    with pytest.raises(FileNotFoundError):
        ph.read(1)
    # ファイルの大きさが量の数と合わなければ止まる
    state = root / "00000000" / "state.0000.f32"
    state.write_bytes(state.read_bytes()[:-4])
    with pytest.raises(ValueError):
        ph.read(0)
    # 知らない版は拒否する
    text = (root / "format.toml").read_text()
    (root / "format.toml").write_text(text.replace("version = 1", "version = 3"))
    with pytest.raises(ValueError):
        pyR2D2.Photosphere(root)
    # 版 1 なのに 6 量でなければ拒否する
    (root / "format.toml").write_text(text.replace('"velocity_horizontal_y"]', '"velocity_horizontal_y", "magnetic_vertical"]'))
    with pytest.raises(ValueError):
        pyR2D2.Photosphere(root)
