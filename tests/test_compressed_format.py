"""R2D2plus の新出力形式 (format = "compressed", R2D2plus DEC-576) の読み込みテスト。

- r2d2plus-z コンテナを Python で書いて read_r2d2plus_z で読み戻す。
- 同じ状態から「従来形式」と「新形式」の小さな合成ランを作り、公開 API
  (d.qf, d.qx, d.qz, d.qm, d.qr, d.qt, d.qs, d.qp, d.qa, d.q2) で読んだ値が
  一致することを確かめる。保存されている変数はビット単位、新形式で計算する
  pr/te/op は独立な倍精度の参照 (pyR2D2.util.eos_switch / eos_table) と
  単精度の丸めの範囲で一致すること。
- 実際の R2D2plus の対照ラン (legacy/compressed) があればビット単位の一致を確かめる
  (無ければ skip)。
"""

import os
from pathlib import Path

import numpy as np
import pytest
from scipy.io import FortranFile

import pyR2D2
from pyR2D2.data_io import compressed as C
from pyR2D2.util import eos_switch, eos_table

numcodecs = pytest.importorskip("numcodecs")

VARS8 = ["ro", "vx", "vy", "vz", "bx", "by", "bz", "se"]
RSTAR = 69598945280.0


# ---------------------------------------------------------------------------
# r2d2plus-z の書き手 (テスト用。R2D2plus compressed_file.cpp と同じ形)
# ---------------------------------------------------------------------------


def _xxh3(raw):
    try:
        import xxhash
    except ImportError:
        return "0" * 16
    return f"{xxhash.xxh3_64_intdigest(raw):016x}"


def write_r2d2plus_z(path, kind, attributes, variables, level=3):
    """variables: list of (name, array, order)。dtype は array の dtype。"""
    blocks, entries, offset = [], [], 0
    for name, array, order in variables:
        array = np.asarray(array)
        dtype = {np.dtype("float32"): "float32", np.dtype("float64"): "float64"}[array.dtype]
        raw = np.ascontiguousarray(array.ravel(order=order)).astype(array.dtype.newbyteorder("<")).tobytes()
        es = array.dtype.itemsize
        shuffled = np.frombuffer(raw, np.uint8).reshape(-1, es).T.tobytes()
        stored = bytes(numcodecs.Zstd(level=level).encode(shuffled))
        blocks.append(stored)
        entries.append(
            "\n[[variables]]\n"
            f'name = "{name}"\ndtype = "{dtype}"\norder = "{order}"\n'
            f"shape = [{', '.join(str(s) for s in array.shape)}]\n"
            'codec = "shuffle+zstd"\n'
            f"offset = {offset}\nstored_bytes = {len(stored)}\nraw_bytes = {len(raw)}\n"
            f'xxh3 = "{_xxh3(raw)}"\n'
        )
        offset += len(stored)
    header = f'format = "r2d2plus-z"\nversion = 1\nkind = "{kind}"\n\n[attributes]\n'
    header += "".join(f"{k} = {v}\n" for k, v in attributes.items())
    header += "".join(entries)
    header = header.encode()
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "wb") as f:
        f.write(C.R2D2PLUS_Z_MAGIC)
        f.write(f"{len(header):016x}\n".encode())
        f.write(header)
        for b in blocks:
            f.write(b)


def test_read_r2d2plus_z_roundtrip(tmp_path):
    rng = np.random.default_rng(0)
    a = rng.standard_normal((3, 4, 5)).astype(np.float32)
    b = rng.standard_normal((7,)).astype(np.float64)
    c = rng.standard_normal((2, 3)).astype(np.float32)
    path = tmp_path / "x.z"
    write_r2d2plus_z(path, "remap_qq", {"nd": 3, "rank": 1},
                     [("a", a, "F"), ("b", b, "C"), ("c", c, "C")])

    header, arrays = C.read_r2d2plus_z(path)
    assert header["kind"] == "remap_qq"
    assert header["attributes"] == {"nd": 3, "rank": 1}
    assert arrays["a"].dtype == np.float32 and arrays["b"].dtype == np.float64
    np.testing.assert_array_equal(arrays["a"], a)
    np.testing.assert_array_equal(arrays["b"], b)
    np.testing.assert_array_equal(arrays["c"], c)

    _, only = C.read_r2d2plus_z(path, names=["c"])
    assert list(only) == ["c"]
    with pytest.raises(KeyError):
        C.read_r2d2plus_z(path, names=["nope"])


def test_read_r2d2plus_z_rejects_other_files(tmp_path):
    path = tmp_path / "raw.bin"
    np.zeros(8, np.float32).tofile(path)
    with pytest.raises(ValueError, match="not a r2d2plus-z"):
        C.read_r2d2plus_z(path)


def test_read_format_toml_absent_means_legacy(tmp_path):
    (tmp_path / "param").mkdir()
    assert C.read_format_toml(tmp_path) is None


def test_read_restart_state_raw_and_compressed(tmp_path):
    nx, ny, nz = 2, 3, 4
    rng = np.random.default_rng(1)
    state = rng.standard_normal(9 * nx * ny * nz)
    meta = (
        '[restart]\nformat_version = 2\ngeneration = ".generation-abc123"\n'
        f"[grid]\nnx = {2 * nx}\nny = {ny}\nnz = {nz}\n"
        "[mpi]\nsize_x = 2\nsize_y = 1\nsize_z = 1\n"
    )
    for compressed in [False, True]:
        slot = tmp_path / ("z" if compressed else "raw")
        gen = slot / ".generation-abc123"
        gen.mkdir(parents=True)
        (slot / "meta.toml").write_text(meta)
        if compressed:
            write_r2d2plus_z(gen / "state.00000001.bin", "restart_state", {},
                             [("state", state, "C")])
        else:
            state.astype("<f8").tofile(gen / "state.00000001.bin")
        out, m = C.read_restart_state(slot, 1)
        assert list(out) == C.RESTART_COMPONENTS
        assert out["ro"].shape == (nx, ny, nz)
        np.testing.assert_array_equal(out["ps"].ravel(), state[8 * nx * ny * nz:])
        # k (C++ z) が最速
        np.testing.assert_array_equal(out["vx"][0, 0, :], state[nx * ny * nz: nx * ny * nz + nz])
        assert m["restart"]["generation"] == ".generation-abc123"


# ---------------------------------------------------------------------------
# 合成ラン (従来形式と新形式を同じ状態から作る)
# ---------------------------------------------------------------------------


def _write_params(path, entries):
    with open(path, "w") as f:
        for key, value, typ in entries:
            f.write(f"{value} {key} {typ}\n")


def _table():
    log_ro = np.linspace(-14.0, -2.0, 49)
    se = np.linspace(0.8e9, 3.2e9, 41)
    lr, s = np.meshgrid(log_ro, se, indexing="ij")
    # 対数で曲がった (双線形では厳密に表せない) 滑らかな量
    pr = np.exp(1.2 * lr + 0.9e-9 * s + 21.0 + 0.02 * np.sin(lr) * np.cos(1e-9 * s))
    te = np.exp(0.25 * lr + 0.6e-9 * s + 9.0 + 0.03 * np.cos(lr + 1e-9 * s))
    op = np.exp(0.5 * lr + 1.1e-9 * s + 3.0 + 0.01 * lr * lr)
    return {"ro": log_ro, "se": se, "pr": pr, "en": 2 * pr, "te": te, "op": op,
            "dprdro": 3 * pr}


def _write_tbl(path, arrays):
    header = ["r2d2plus-table 1", "kind eos_table_sero",
              "uuid 00000000-0000-0000-0000-000000000000",
              "layout row_major little_endian float64", "binary_offset 1024"]
    for k, a in arrays.items():
        header.append("array " + k + " " + " ".join(str(n) for n in a.shape))
    header.append("end")
    text = ("\n".join(header) + "\n").encode()
    with open(path, "wb") as f:
        f.write(text + b"\n" * (1024 - len(text)))
        for a in arrays.values():
            f.write(np.ascontiguousarray(a, dtype="<f8").tobytes())


class _TableData:
    """util.eos_table / eos_switch 用 (倍精度の独立な参照)。"""

    def __init__(self, table, back):
        self.p = self
        self.log_ro_e = table["ro"]
        self.se_e = table["se"]
        self.dlogro_e = self.log_ro_e[1] - self.log_ro_e[0]
        self.dse_e = self.se_e[1] - self.se_e[0]
        for q in ["pr", "te", "op"]:
            setattr(self, f"log_{q}_e", np.log(table[q] + 1.0e-200))
        for k, v in back.items():
            setattr(self, k, v)


def _background(table, x):
    ix = x.size
    ro0 = np.exp(np.linspace(-6.0, -9.0, ix))
    se0 = np.linspace(1.9e9, 2.1e9, ix)
    ref = _TableData(table, {})
    pr0 = eos_table(ref, ro0, se0, "pr")
    te0 = eos_table(ref, ro0, se0, "te")
    back = {"ro0": ro0, "se0": se0, "pr0": pr0, "te0": te0,
            "dprdro": 1.3 * pr0 / ro0, "dprdse": 0.9e-9 * pr0,
            "dtedro": 0.2 * te0 / ro0, "dtedse": 0.6e-9 * te0}
    return back


def _write_back(path, xg, yg, zg, back, margin):
    ixg = xg.size
    names = ["pr0", "te0", "ro0", "se0", "en0", "op0", "tu0", "dsedr0", "dtedr0",
             "dprdro", "dprdse", "dtedro", "dtedse", "dendro", "dendse", "gx", "cp",
             "fa", "sa", "xi"]
    with open(path, "wb") as f:
        np.zeros(1, "<i4").tofile(f)
        for a in (xg, yg, zg):
            a.astype("<f8").tofile(f)
        for name in names:
            v = np.zeros(ixg)
            if name in back:
                v[margin:ixg - margin] = back[name]
                v[:margin] = back[name][0]
                v[ixg - margin:] = back[name][-1]
            v.astype("<f8").tofile(f)
        np.zeros(1, "<i4").tofile(f)


FORMAT_TOML = """format = "compressed"
version = 1
variables = ["ro", "vx", "vy", "vz", "bx", "by", "bz", "se"]
derived_not_written = ["pr", "te", "op"]
eos_table = "eos_table.tbl"
rstar = 69598945280
compression_level = 3
remap_qq = "remap/qq/{group:05d}/{rank:08d}/qq.z.{nd:08d}.{rank:08d}"
tau = "tau/qq.c.{nd_tau:08d}"
tau_quantities = ["in", "ro", "se", "vx", "vy", "vz", "bx", "by", "bz", "height", "fr"]
tau_height_origin = "rstar"
slice = "slice/{name}.c.{nd_tau:08d}.{n:08d}"
remap_2d = "remap/qq/qq.c.{nd:08d}"
remap_2d_extra = ["tu", "fr"]
prev = "prev/{group:05d}/{rank:08d}/qq.z.{nd:08d}.{index:08d}.{rank:08d}"
aftr = "aftr/{group:05d}/{rank:08d}/qq.z.{nd:08d}.{index:08d}.{rank:08d}"
"""


def _perturb(rng, shape, ro0, se0):
    """半分ほどの点が ct0 = 3e-3 を超える摂動 (単精度で表せる値)。"""
    amp = np.where(rng.random(shape) < 0.5, 1e-3, 1e-2)
    ro = (ro0 * amp * rng.uniform(-1, 1, shape)).astype(np.float32)
    se = (se0 * amp * rng.uniform(-1, 1, shape)).astype(np.float32)
    return ro, se


def _make_3d_runs(root):
    """同じ状態から従来形式と新形式の 3D ランを作る。戻り値は (legacy, compressed, truth)。"""
    nx, ny, nz, ix0, jx0, kx0, margin = 8, 4, 4, 1, 2, 1, 2
    npe = ix0 * jx0 * kx0
    ix, jx, kx = ix0 * nx, jx0 * ny, kx0 * nz
    ixg, jxg, kxg = ix + 2 * margin, jx + 2 * margin, kx + 2 * margin
    ixr, jxr = 2, 1
    rng = np.random.default_rng(2)
    table = _table()
    dx = 4.8e7
    xg = RSTAR - 3e8 + dx * (np.arange(ixg) - margin + 0.5)
    yg = dx * (np.arange(jxg) - margin + 0.5)
    zg = dx * (np.arange(kxg) - margin + 0.5)
    x = xg[margin:-margin]
    back = _background(table, x)
    ref = _TableData(table, back)

    state = {}
    ro0 = back["ro0"][:, None, None]
    se0 = back["se0"][:, None, None]
    state["ro"], state["se"] = _perturb(rng, (ix, jx, kx), ro0, se0)
    for k in ["vx", "vy", "vz", "bx", "by", "bz", "ph"]:
        state[k] = rng.standard_normal((ix, jx, kx)).astype(np.float32) * 1e5
    ro64, se64 = state["ro"].astype(np.float64), state["se"].astype(np.float64)
    truth = {k: eos_switch(ref, ro64, se64, k).astype(np.float32) for k in ["pr", "te", "op"]}

    runs = {}
    for fmt in ["legacy", "compressed"]:
        d = root / fmt / "data"
        (d / "param").mkdir(parents=True)
        _write_params(d / "param" / "params.dac", [
            ("xdcheck", 2, "i"), ("ydcheck", 2, "i"), ("zdcheck", 2, "i"),
            ("margin", margin, "i"), ("nx", nx, "i"), ("ny", ny, "i"), ("nz", nz, "i"),
            ("npe", npe, "i"), ("ix0", ix0, "i"), ("jx0", jx0, "i"), ("kx0", kx0, "i"),
            ("mtype", 9, "i"), ("swap", 0, "i"), ("ixr", ixr, "i"), ("jxr", jxr, "i"),
            ("m_in", 13, "i"), ("m_tu", 3, "i"), ("rstar", RSTAR, "d"),
            ("geometry", "Cartesian", "c"), ("ib_rte_bot", 0, "i"),
        ])
        (d / "param" / "nd.dac").write_text("1 1\n")
        ff = FortranFile(d / "param" / "xyz.dac", "w")
        xyz = np.array([[0, j, 0] for j in range(jx0)], dtype=np.int32)
        ff.write_record(xyz.ravel(order="F"))
        ff.close()
        _write_back(d / "param" / "back.dac", xg, yg, zg, back, margin)
        (d / "remap" / "vl").mkdir(parents=True)
        (d / "remap" / "vl" / "c.dac").write_text("0 0 0\n")
        info = {
            "iss": [1, nx // 2 + 1], "iee": [nx // 2, nx], "jss": [1, 1], "jee": [jx, jx],
            "iixl": [nx // 2, nx // 2], "jjxl": [jx, jx], "np_ijr": [0, 1],
            "ir": [1, 2], "jr": [1, 1],
            "i2ir": [1] * (margin + nx // 2) + [2] * (margin + nx // 2),
            "j2jr": [1] * jxg,
        }
        with open(d / "remap" / "remap_info.dac", "wb") as f:
            for key in ["iss", "iee", "jss", "jee", "iixl", "jjxl", "np_ijr", "ir", "jr",
                        "i2ir", "j2jr"]:
                np.asarray(info[key], "<i4").tofile(f)

        n = 1
        for rank in range(ixr):
            i0, i1 = info["iss"][rank] - 1, info["iee"][rank]
            block = {k: np.asfortranarray(v[i0:i1]) for k, v in {**state, **truth}.items()}
            sub = d / "remap" / "qq" / "00000" / f"{rank:08d}"
            sub.mkdir(parents=True)
            if fmt == "legacy":
                with open(sub / f"qq.dac.{n:08d}.{rank:08d}", "wb") as f:
                    for k in VARS8 + ["ph", "pr", "te", "op"]:
                        block[k].ravel(order="F").astype("<f4").tofile(f)
            else:
                write_r2d2plus_z(sub / f"qq.z.{n:08d}.{rank:08d}", "remap_qq",
                                 {"nd": n, "rank": rank},
                                 [(k, block[k], "F") for k in VARS8])

        # tau: 面の値 (全量)。1 本だけ「面が見つからない」柱 (全部 0) を入れる
        m_tu = 3
        tau = {}
        itau = [nx - 3, nx - 2, nx - 1]
        for mt, i in enumerate(itau):
            trng = np.random.default_rng(10 + mt)  # 両形式で同じ値にする
            s = {k: state[k][i].astype(np.float32) for k in state}
            ro_t = (s["ro"].astype(np.float64) + back["ro0"][i]).astype(np.float32)
            se_t = (s["se"].astype(np.float64) + back["se0"][i]).astype(np.float32)
            r = x[i] + 1.234e5 * trng.random((jx, kx))
            fr = trng.random((jx, kx)).astype(np.float32)
            tau[mt] = dict(ro=ro_t, se=se_t, vx=s["vx"], vy=s["vy"], vz=s["vz"],
                           bx=s["bx"], by=s["by"], bz=s["bz"], fr=fr, r=r,
                           rt=(trng.random((jx, kx)) if mt == 0 else np.zeros((jx, kx))).astype(np.float32))
            for key in ["ro", "se"]:
                tau[mt][key][0, 0] = 0.0
            for key in ["pr", "te"]:
                with np.errstate(all="ignore"):
                    v = eos_table(ref, tau[mt]["ro"].astype(np.float64),
                                  tau[mt]["se"].astype(np.float64), key).astype(np.float32)
                v[0, 0] = 0.0
                tau[mt][key] = v
        (d / "tau").mkdir()
        if fmt == "legacy":
            names = ["rt", "ro", "se", "pr", "te", "vx", "vy", "vz", "bx", "by", "bz", "he", "fr"]
            buf = np.zeros((m_tu, len(names), jx, kx), np.float32)
            for mt in range(m_tu):
                for m, k in enumerate(names):
                    buf[mt, m] = tau[mt]["r"].astype(np.float32) if k == "he" else tau[mt][k]
            buf.ravel(order="F").astype("<f4").tofile(d / "tau" / f"qq.dac.{n:08d}")
        else:
            names = ["rt", "ro", "se", "vx", "vy", "vz", "bx", "by", "bz", "height", "fr"]
            buf = np.zeros((m_tu, len(names), jx, kx), np.float32)
            for mt in range(m_tu):
                for m, k in enumerate(names):
                    buf[mt, m] = (tau[mt]["r"] - RSTAR).astype(np.float32) if k == "height" else tau[mt][k]
            buf.ravel(order="F").astype("<f4").tofile(d / "tau" / f"qq.c.{n:08d}")

        # slice: 各方向 1 枚
        (d / "slice").mkdir()
        (d / "slice" / "params.dac").write_text(
            "       1 nx_slice i\n       1 ny_slice i\n       1 nz_slice i\n")
        i_s, j_s, k_s = nx - 2, 1, 2
        np.array([x[i_s], yg[margin + j_s], zg[margin + k_s]], "<f8").tofile(d / "slice" / "slice.dac")
        planes = {
            "x": ({k: v[i_s] for k, v in state.items()}, i_s),
            "y": ({k: v[:, j_s, :] for k, v in state.items()}, None),
            "z": ({k: v[:, :, k_s] for k, v in state.items()}, None),
        }
        for direc, (s, i_fixed) in planes.items():
            if i_fixed is None:
                s_ro, s_se = s["ro"].astype(np.float64), s["se"].astype(np.float64)
                prte = {k: eos_switch(ref, s_ro, s_se, k).astype(np.float32) for k in ["pr", "te"]}
            else:
                # 動径の断面: 鉛直位置は 1 点。背景をその高さに固定して参照を作る
                one = _TableData(table, {k: np.full(1, v[i_fixed]) for k, v in back.items()})
                s_ro = s["ro"].astype(np.float64)[None]
                s_se = s["se"].astype(np.float64)[None]
                prte = {k: eos_switch(one, s_ro, s_se, k)[0].astype(np.float32) for k in ["pr", "te"]}
            keys = VARS8 + ["ph", "pr", "te"] if fmt == "legacy" else VARS8
            allv = {**s, **prte}
            buf = np.stack([allv[k] for k in keys], axis=-1)
            name = f"qq{direc}.dac.{n:08d}.{1:08d}" if fmt == "legacy" else f"qq{direc}.c.{n:08d}.{1:08d}"
            buf.ravel(order="F").astype("<f4").tofile(d / "slice" / name)

        # prev / aftr (従来は倍精度、新形式は単精度)
        for kind in ["prev", "aftr"]:
            for rank in range(npe):
                jb = rank
                blk = {k: state[k][:, jb * ny:(jb + 1) * ny, :] for k in VARS8}
                sub = d / kind / "00000" / f"{rank:08d}"
                sub.mkdir(parents=True)
                if fmt == "legacy":
                    buf = np.stack([blk[k].astype(np.float64) for k in VARS8], axis=-1)
                    buf.ravel(order="F").astype("<f8").tofile(sub / f"qq.dac.{n:08d}.{0:08d}.{rank:08d}")
                else:
                    write_r2d2plus_z(sub / f"qq.z.{n:08d}.{0:08d}.{rank:08d}", f"raw_{kind}",
                                     {"index": 0, "nd": n, "rank": rank},
                                     [(k, np.asfortranarray(blk[k]), "F") for k in VARS8])

        if fmt == "compressed":
            (d / "param" / "format.toml").write_text(FORMAT_TOML)
            _write_tbl(d / "param" / "eos_table.tbl", table)
        runs[fmt] = d
    return runs["legacy"], runs["compressed"], dict(state=state, truth=truth, tau=tau)


@pytest.fixture(scope="module")
def runs3d(tmp_path_factory):
    legacy, comp, truth = _make_3d_runs(tmp_path_factory.mktemp("runs3d"))
    return pyR2D2.Data(legacy), pyR2D2.Data(comp), truth


def _bitwise(a, b):
    a, b = np.asarray(a), np.asarray(b)
    assert a.shape == b.shape
    assert a.dtype == b.dtype
    np.testing.assert_array_equal(
        np.ascontiguousarray(a).view(np.uint8), np.ascontiguousarray(b).view(np.uint8)
    )


def _close_f32(computed, reference, scale):
    """単精度の丸め (と線形 EOS の単精度演算) の範囲で一致すること。"""
    computed = np.asarray(computed, np.float64)
    reference = np.asarray(reference, np.float64)
    assert computed.shape == reference.shape
    np.testing.assert_array_less(np.abs(computed - reference), 1e-6 * np.abs(scale) + 1e-30)


def test_detection(runs3d):
    L, Cd, _ = runs3d
    assert L.p.output_format == "legacy" and L.p.format_info is None
    assert Cd.p.output_format == "compressed"
    assert Cd.p.format_info["variables"] == VARS8
    # EOS 表は param/eos_table.tbl から自動で読まれる
    assert Cd.p.eos_table_path.endswith("eos_table.tbl")
    assert Cd.log_ro_e.size == 49


def test_fulldata(runs3d):
    L, Cd, t = runs3d
    L.qf.read(1)
    Cd.qf.read(1)
    for k in VARS8:
        _bitwise(Cd.qf.__dict__[k], L.qf.__dict__[k])
    # 従来形式のファイルに入れた pr/te/op (倍精度の独立な参照) と一致
    scales = {
        "pr": np.broadcast_to(Cd.pr0[:, None, None], L.qf.pr.shape),
        "te": np.broadcast_to(Cd.te0[:, None, None], L.qf.te.shape),
        "op": L.qf.op,
    }
    for k in ["pr", "te", "op"]:
        assert Cd.qf.__dict__[k].dtype == np.float32
        _close_f32(Cd.qf.__dict__[k], L.qf.__dict__[k], scales[k])
    # 表と線形の両方の枝を踏んでいること
    ct = np.maximum(np.abs(t["state"]["ro"]) / Cd.ro0[:, None, None],
                    np.abs(t["state"]["se"]) / Cd.se0[:, None, None])
    assert 0.2 < (ct >= 3e-3).mean() < 0.8
    # ph は書かれていないので "all" では None
    assert Cd.qf.ph is None


def test_fulldata_keys_and_missing_ph(runs3d):
    L, Cd, _ = runs3d
    Cd.qf.read(1, keys=["te"])
    L.qf.read(1, keys=["te"])
    _close_f32(Cd.qf.te, L.qf.te, np.broadcast_to(Cd.te0[:, None, None], L.qf.te.shape))
    with pytest.raises(ValueError, match="ph"):
        Cd.qf.read(1, keys=["ph"])


def test_plane_and_region_readers(runs3d):
    L, Cd, _ = runs3d
    xs = L.x[5]
    L.qx.read(xs, 1)
    Cd.qx.read(xs, 1)
    L.qz.read(L.z[2], 1)
    Cd.qz.read(L.z[2], 1)
    L.qm.read(1, 1)
    Cd.qm.read(1, 1)
    kw = dict(x0=L.x[2], x1=L.x[6], y0=L.y[1], y1=L.y[6], z0=L.z[1], z1=L.z[2])
    L.qr.read(1, ["ro", "vz", "pr", "te", "op"], **kw)
    Cd.qr.read(1, ["ro", "vz", "pr", "te", "op"], **kw)
    for r_l, r_c in [(L.qx, Cd.qx), (L.qz, Cd.qz), (L.qm, Cd.qm), (L.qr, Cd.qr)]:
        for k in ["ro", "vz", "se"] if r_l is not L.qr else ["ro", "vz"]:
            _bitwise(r_c.__dict__[k], r_l.__dict__[k])
        for k in ["pr", "te", "op"]:
            ref = r_l.__dict__[k]
            scale = {"pr": Cd.pr0.max(), "te": Cd.te0.max(), "op": np.abs(ref)}[k]
            _close_f32(r_c.__dict__[k], ref, scale)
    # 部分読みは全体読みの切り出しとビット単位で同じ
    Cd.qf.read(1)
    _bitwise(Cd.qx.te, Cd.qf.te[5])
    _bitwise(Cd.qz.pr, Cd.qf.pr[:, :, 2])
    _bitwise(Cd.qr.op, Cd.qf.op[2:7, 1:7, 1:3])


def test_optical_depth(runs3d):
    L, Cd, t = runs3d
    L.qt.read(1)
    Cd.qt.read(1)
    for tau in ["", "01", "001"]:
        for k in ["rt", "ro", "se", "vx", "vy", "vz", "bx", "by", "bz", "fr"]:
            _bitwise(Cd.qt.__dict__[k + tau], L.qt.__dict__[k + tau])
        for k in ["pr", "te"]:
            ref = L.qt.__dict__[k + tau]
            _close_f32(Cd.qt.__dict__[k + tau], ref, np.abs(ref))
            assert Cd.qt.__dict__[k + tau][0, 0] == 0.0  # 面が無い柱は 0
        # he は中心からの半径 (倍精度)、height は r - rstar (単精度のまま)
        he_c = Cd.qt.__dict__["he" + tau]
        assert he_c.dtype == np.float64
        assert np.abs(he_c - L.qt.__dict__["he" + tau]).max() <= 4096.0
        np.testing.assert_allclose(Cd.qt.__dict__["height" + tau], he_c - RSTAR, rtol=0, atol=1e-6)


def test_slice(runs3d):
    L, Cd, _ = runs3d
    for direc in ["x", "y", "z"]:
        L.qs.read(0, direc, 1)
        Cd.qs.read(0, direc, 1)
        for k in VARS8:
            _bitwise(Cd.qs.__dict__[k], L.qs.__dict__[k])
        for k in ["pr", "te"]:
            ref = L.qs.__dict__[k]
            _close_f32(Cd.qs.__dict__[k], ref, np.maximum(np.abs(ref), 1e-3 * np.abs(ref).max()))


def test_prev_aftr(runs3d):
    L, Cd, _ = runs3d
    for a, b in [(L.qp, Cd.qp), (L.qa, Cd.qa)]:
        a.read(1, 0)
        b.read(1, 0)
        for k in VARS8:
            assert b.__dict__[k].dtype == np.float32
            assert b.__dict__[k].shape == a.__dict__[k].shape
            # 従来の倍精度は単精度の値をそのまま広げたもの
            _bitwise(b.__dict__[k], a.__dict__[k].astype(np.float32))


def test_missing_eos_table_gives_clear_error(tmp_path):
    _, comp, _ = _make_3d_runs(tmp_path)
    (comp / "param" / "eos_table.tbl").unlink()
    d = pyR2D2.Data(comp)
    d.qf.read(1, keys=["ro", "se"])  # 保存量だけなら読める
    with pytest.raises(FileNotFoundError, match="eos_table"):
        d.qf.read(1, keys=["pr"])


# ---------------------------------------------------------------------------
# 2D (ny = 1 相当、Fortran の zdcheck = 1)
# ---------------------------------------------------------------------------


def _make_2d_runs(root):
    nx, ny, margin = 8, 6, 2
    ix, jx = nx, ny
    ixg, jxg = ix + 2 * margin, jx + 2 * margin
    rng = np.random.default_rng(3)
    table = _table()
    dx = 4.8e7
    xg = RSTAR - 3e8 + dx * (np.arange(ixg) - margin + 0.5)
    yg = dx * (np.arange(jxg) - margin + 0.5)
    zg = np.zeros(1)
    back = _background(table, xg[margin:-margin])
    ref = _TableData(table, back)
    state = {}
    state["ro"], state["se"] = _perturb(rng, (ix, jx), back["ro0"][:, None], back["se0"][:, None])
    for k in ["vx", "vy", "vz", "bx", "by", "bz"]:
        state[k] = rng.standard_normal((ix, jx)).astype(np.float32)
    state["tu"] = rng.random((ix, jx)).astype(np.float32)
    state["fr"] = rng.random((ix, jx)).astype(np.float32)
    rt = state["ro"].astype(np.float64) + back["ro0"][:, None]
    st = state["se"].astype(np.float64) + back["se0"][:, None]
    # 従来形式の 2D は切り替え無しの表の値 (R2D2plus remap2d_writer.cpp)
    state["pr"] = (eos_table(ref, rt, st, "pr") - back["pr0"][:, None]).astype(np.float32)
    state["te"] = (eos_table(ref, rt, st, "te") - back["te0"][:, None]).astype(np.float32)
    state["op"] = eos_table(ref, rt, st, "op").astype(np.float32)
    runs = {}
    for fmt in ["legacy", "compressed"]:
        d = root / fmt / "data"
        (d / "param").mkdir(parents=True)
        _write_params(d / "param" / "params.dac", [
            ("xdcheck", 2, "i"), ("ydcheck", 2, "i"), ("zdcheck", 1, "i"),
            ("margin", margin, "i"), ("nx", nx, "i"), ("ny", ny, "i"), ("nz", 1, "i"),
            ("npe", 1, "i"), ("ix0", 1, "i"), ("jx0", 1, "i"), ("kx0", 1, "i"),
            ("mtype", 9, "i"), ("swap", 0, "i"), ("rstar", RSTAR, "d"),
            ("geometry", "Cartesian", "c"),
        ])
        (d / "param" / "nd.dac").write_text("1 1\n")
        ff = FortranFile(d / "param" / "xyz.dac", "w")
        ff.write_record(np.zeros(3, np.int32))
        ff.close()
        _write_back(d / "param" / "back.dac", xg, yg, zg, back, margin)
        (d / "remap" / "qq").mkdir(parents=True)
        if fmt == "legacy":
            names = VARS8 + ["pr", "te", "op", "tu", "fr"]
            buf = np.zeros((14, ix, jx), np.float32)
            for m, k in enumerate(names):
                buf[m] = state[k]
            buf.ravel(order="F").astype("<f4").tofile(d / "remap" / "qq" / "qq.dac.00000001")
        else:
            names = VARS8 + ["tu", "fr"]
            buf = np.stack([state[k] for k in names])
            buf.ravel(order="F").astype("<f4").tofile(d / "remap" / "qq" / "qq.c.00000001")
            (d / "param" / "format.toml").write_text(FORMAT_TOML)
            _write_tbl(d / "param" / "eos_table.tbl", table)
        runs[fmt] = d
    return runs["legacy"], runs["compressed"], state


def test_two_dimension(tmp_path):
    legacy, comp, state = _make_2d_runs(tmp_path)
    L, Cd = pyR2D2.Data(legacy), pyR2D2.Data(comp)
    L.q2.read(1)
    Cd.q2.read(1)
    for k in VARS8:
        _bitwise(Cd.q2.__dict__[k], L.q2.__dict__[k])
    # 従来の読み手は pr まで (添字 8) しか埋めないので、残りはファイルの値と比べる
    _close_f32(Cd.q2.pr, L.q2.pr, np.broadcast_to(Cd.pr0[:, None], L.q2.pr.shape))
    _close_f32(Cd.q2.te, state["te"], np.broadcast_to(Cd.te0[:, None], state["te"].shape))
    _close_f32(Cd.q2.op, state["op"], state["op"])
    _bitwise(Cd.q2.tu, state["tu"])
    _bitwise(Cd.q2.fr, state["fr"])
    for k in Cd.q2.value_keys:
        assert Cd.q2.__dict__[k].shape == (L.ix, L.jx)


# ---------------------------------------------------------------------------
# 実ランの対照 (R2D2plus の legacy / compressed を同じ設定で走らせたもの)
# ---------------------------------------------------------------------------

FIXTURE = Path(os.environ.get("PYR2D2_IO_FIXTURE", "/scr/a000/c0234hotta/io_fixture_20260927"))


@pytest.mark.skipif(
    not (FIXTURE / "comp_3d" / "data" / "param" / "format.toml").exists(),
    reason="R2D2plus I/O fixture not available",
)
def test_r2d2plus_fixture_matches_legacy():
    L = pyR2D2.Data(FIXTURE / "legacy_3d" / "data")
    Cd = pyR2D2.Data(FIXTURE / "comp_3d" / "data")
    assert Cd.p.output_format == "compressed"
    n = 6
    L.qf.read(n)
    Cd.qf.read(n)
    # 保存量も pr/te/op もビット単位で一致 (R2D2plus と同じ式・丸め)
    for k in VARS8 + ["pr", "te", "op"]:
        _bitwise(Cd.qf.__dict__[k], L.qf.__dict__[k])
    L.qt.read(n)
    Cd.qt.read(n)
    for tau in ["", "01", "001"]:
        for k in ["rt", "ro", "se", "pr", "te", "vx", "vy", "vz", "bx", "by", "bz", "fr"]:
            if k == "rt" and tau:
                continue
            _bitwise(Cd.qt.__dict__[k + tau], L.qt.__dict__[k + tau])
        assert np.abs(Cd.qt.__dict__["he" + tau] - L.qt.__dict__["he" + tau]).max() <= 4096.0
        _bitwise(Cd.qt.__dict__["he" + tau].astype(np.float32), L.qt.__dict__["he" + tau])
    for direc in ["x", "y", "z"]:
        L.qs.read(0, direc, n)
        Cd.qs.read(0, direc, n)
        for k in VARS8:
            _bitwise(Cd.qs.__dict__[k], L.qs.__dict__[k])
        # slice の pr/te は R2D2plus が倍精度の状態から作ったもの。単精度の入力からの
        # 再計算なので、背景量に対して 1e-7 以内 (実測 2e-8)
        for k, X0 in [("pr", Cd.pr0), ("te", Cd.te0)]:
            if direc == "x":
                X0 = X0[np.argmin(np.abs(Cd.x - Cd.x_slice[0]))]
            else:
                X0 = X0[:, None]
            diff = np.abs(Cd.qs.__dict__[k].astype(float) - L.qs.__dict__[k])
            assert (diff / np.abs(X0)).max() < 1e-7
    L.qp.read(1, 0)
    Cd.qp.read(1, 0)
    for k in VARS8:
        _bitwise(Cd.qp.__dict__[k], L.qp.__dict__[k].astype(np.float32))
    state, meta = C.read_restart_state(FIXTURE / "comp_3d" / "data" / "restart" / "e", 0)
    assert state["ro"].shape == (64, 64, 128)
    assert meta["restart"]["nd"] == 7


# ---------------------------------------------------------------------------
# 初期状態ファイルの書き出し (write_initial_state)
# ---------------------------------------------------------------------------


def test_write_initial_state_roundtrip(tmp_path):
    pytest.importorskip("xxhash")
    rng = np.random.default_rng(4)
    shape = (6, 5, 4)  # (ix, jx, kx) = (鉛直, 第1水平, 第2水平)
    fields = {k: rng.standard_normal(shape) for k in VARS8}
    fields["vx"] = np.asfortranarray(fields["vx"])  # 並びの違う入力でもよい
    fields["se"] = fields["se"].astype(np.float32)  # float64 に広げて書く
    attrs = {"source": 'run "d001"', "nd": 120, "time": 3600.0, "restart": False}
    path = pyR2D2.write_initial_state(tmp_path / "init.z", fields, attributes=attrs)

    header, arrays = C.read_r2d2plus_z(path, verify=True)  # XXH3 も照合する
    assert header["kind"] == "initial_state"
    assert header["attributes"] == attrs
    assert [v["name"] for v in header["variables"]] == VARS8
    for v in header["variables"]:
        assert v["dtype"] == "float64" and v["order"] == "F" and v["shape"] == list(shape)
    for k in VARS8:
        _bitwise(arrays[k], np.asarray(fields[k], dtype=np.float64))
    # 鉛直 (第 0 軸) が最速で並んでいること: 生バイトを直接見る
    raw = C._unshuffle(C._zstd_decompress(
        open(path, "rb").read()[header["data_offset"]:][:header["variables"][0]["stored_bytes"]],
        header["variables"][0]["raw_bytes"]), "<f8")
    np.testing.assert_array_equal(raw[: shape[0]], fields["ro"][:, 0, 0])


def test_write_initial_state_rejects_bad_input(tmp_path):
    pytest.importorskip("xxhash")
    good = np.zeros((4, 3, 2))
    with pytest.raises(ValueError, match="unknown"):
        pyR2D2.write_initial_state(tmp_path / "a.z", {"pr": good})
    with pytest.raises(ValueError, match="shape"):
        pyR2D2.write_initial_state(tmp_path / "a.z", {"ro": good, "se": np.zeros((4, 3, 3))})
    with pytest.raises(ValueError, match="3D"):
        pyR2D2.write_initial_state(tmp_path / "a.z", {"ro": np.zeros((4, 3))})
    with pytest.warns(UserWarning, match="perturbation"):
        pyR2D2.write_initial_state(tmp_path / "a.z", {"ro": good + 1e-7})


def test_write_initial_state_needs_xxhash(tmp_path, monkeypatch):
    import builtins

    original_import = builtins.__import__

    def no_xxhash(name, *args, **kwargs):
        if name == "xxhash":
            raise ImportError("blocked by test")
        return original_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", no_xxhash)
    with pytest.raises(ImportError, match="xxhash"):
        pyR2D2.write_initial_state(tmp_path / "a.z", {"ro": np.zeros((2, 2, 2))})


@pytest.mark.skipif(
    not (FIXTURE / "comp_3d" / "data" / "param" / "format.toml").exists(),
    reason="R2D2plus I/O fixture not available",
)
def test_fixture_checksums_and_initial_state_from_output(tmp_path):
    pytest.importorskip("xxhash")
    # xxhash の XXH3-64 が R2D2plus (XXH3_64bits) と同じであること
    C.read_r2d2plus_z(
        FIXTURE / "comp_3d" / "data" / "remap" / "qq" / "00000" / "00000002"
        / "qq.z.00000004.00000002", verify=True)
    # 前の出力を初期条件にする (従来形式のランから)
    d = pyR2D2.Data(FIXTURE / "legacy_3d" / "data")
    d.qf.read(6, keys=VARS8)
    fields = {k: d.qf.__dict__[k] for k in VARS8}
    path = pyR2D2.write_initial_state(tmp_path / "init.z", fields, attributes={"nd": 6})
    _, arrays = C.read_r2d2plus_z(path, verify=True)
    for k in VARS8:
        assert arrays[k].shape == (d.ix, d.jx, d.kx)
        _bitwise(arrays[k], d.qf.__dict__[k].astype(np.float64))
