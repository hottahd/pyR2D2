"""C++球面→直交変換の精度、dtype、および入力条件を確認する。"""

import numpy as np
import pytest

from pyR2D2.cpp_util import spherical2cartesian


@pytest.fixture
def spherical_grid():
    """一様な球面座標を返す。"""
    rr = np.linspace(0.0, 5.0, 6)
    th = np.linspace(0.0, np.pi, 9)
    ph = np.linspace(-np.pi, np.pi, 9)
    return rr, th, ph


def linear_radial_colatitude_field(rr, th, ph):
    """三線形補間で再現でき、φに依存しない試験場を返す。"""
    return (
        3.0
        + 2.0 * rr[:, None, None]
        - 0.5 * th[None, :, None]
        + np.zeros_like(ph)[None, None, :]
    )


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
def test_matches_analytic_field_and_preserves_dtype(spherical_grid, dtype):
    """解析場を再現し、入力と同じdtypeのC連続配列を返すことを確認する。"""
    rr, th, ph = spherical_grid
    qqs = np.ascontiguousarray(
        linear_radial_colatitude_field(rr, th, ph),
        dtype=dtype,
    )

    result, xc, yc, zc = spherical2cartesian(
        rr, th, ph, qqs, ph.size, ph.size, ph.size
    )

    xx = xc[:, None, None]
    yy = yc[None, :, None]
    zz = zc[None, None, :]
    radius = np.sqrt(xx**2 + yy**2 + zz**2)
    colatitude = np.arctan2(np.sqrt(xx**2 + yy**2), zz)
    inside = (
        (radius > rr[0])
        & (radius < rr[-2])
        & (colatitude > th[0])
        & (colatitude < th[-2])
    )
    expected = np.zeros(result.shape, dtype=np.float64)
    expected[inside] = (3.0 + 2.0 * radius - 0.5 * colatitude)[inside]

    assert result.dtype == dtype
    assert result.flags.c_contiguous
    assert xc.dtype == yc.dtype == zc.dtype == np.float64
    tolerance = 2.0e-6 if dtype == np.float32 else 1.0e-12
    np.testing.assert_allclose(result, expected, rtol=0.0, atol=tolerance)


def test_phi_period_does_not_depend_on_output_z_size(spherical_grid):
    """φ周期折返しが無関係な出力z格子数ではなく球面φ格子に基づくことを確認する。"""
    rr, th, ph = spherical_grid
    qqs = np.ascontiguousarray(
        np.broadcast_to(
            np.cos(ph)[None, None, :],
            (rr.size, th.size, ph.size),
        ),
        dtype=np.float32,
    )

    result, xc, yc, zc = spherical2cartesian(rr, th, ph, qqs, 9, 9, 7)

    ix = np.flatnonzero(xc == -2.5).item()
    iy = np.flatnonzero(yc == 0.0).item()
    iz = np.argmin(np.abs(zc - 5.0 / 3.0))
    np.testing.assert_allclose(result[ix, iy, iz], -1.0, atol=1.0e-6)


def test_accepts_float32_coordinates(spherical_grid):
    """float32座標を小さなfloat64配列へ変換して利用できることを確認する。"""
    rr, th, ph = (coordinate.astype(np.float32) for coordinate in spherical_grid)
    qqs = np.ones((rr.size, th.size, ph.size), dtype=np.float32)

    result, _, _, _ = spherical2cartesian(rr, th, ph, qqs, 5, 5, 5)

    assert result.dtype == np.float32


@pytest.mark.parametrize("invalid_dtype", [np.int32, np.int64, np.float16])
def test_rejects_unsupported_input_dtype(spherical_grid, invalid_dtype):
    """暗黙の巨大変換を避けるため、float32/float64以外を拒否する。"""
    rr, th, ph = spherical_grid
    qqs = np.ones((rr.size, th.size, ph.size), dtype=invalid_dtype)

    with pytest.raises(TypeError, match="float32 or float64"):
        spherical2cartesian(rr, th, ph, qqs, 5, 5, 5)


def test_accepts_f_contiguous_input(spherical_grid):
    """F連続入力をコピーせず読み、C連続入力と同じ結果を返す。"""
    rr, th, ph = spherical_grid
    qqs_c = np.ascontiguousarray(
        linear_radial_colatitude_field(rr, th, ph), dtype=np.float32
    )
    qqs_f = np.asfortranarray(qqs_c)

    expected, _, _, _ = spherical2cartesian(rr, th, ph, qqs_c, 5, 5, 5)
    result, _, _, _ = spherical2cartesian(rr, th, ph, qqs_f, 5, 5, 5)

    assert qqs_f.flags.f_contiguous
    np.testing.assert_array_equal(result, expected)


def test_rejects_noncontiguous_input(spherical_grid):
    """C連続でもF連続でもない飛び飛びviewを拒否する。"""
    rr, th, ph = spherical_grid
    storage = np.ones((rr.size, th.size, 2 * ph.size), dtype=np.float32)
    qqs = storage[:, :, ::2]

    assert not qqs.flags.c_contiguous
    assert not qqs.flags.f_contiguous
    with pytest.raises(ValueError, match="C- or F-contiguous"):
        spherical2cartesian(rr, th, ph, qqs, 5, 5, 5)


def test_accepts_native_float_dtype_with_metadata(spherical_grid):
    """dtype metadataを持つR2D2配列もコピーせず受け付ける。"""
    rr, th, ph = spherical_grid
    dtype = np.dtype(np.float32, metadata={"source": "r2d2"})
    qqs = np.ones((rr.size, th.size, ph.size), dtype=dtype)

    result, _, _, _ = spherical2cartesian(rr, th, ph, qqs, 5, 5, 5)

    assert result.dtype == np.float32


def test_rejects_source_coordinate_length_mismatch(spherical_grid):
    """球面座標長とqqsの対応軸長が異なる入力を拒否する。"""
    rr, th, ph = spherical_grid
    qqs = np.ones((rr.size, th.size, ph.size), dtype=np.float32)

    with pytest.raises(ValueError, match="spherical coordinate lengths"):
        spherical2cartesian(rr[:-1], th, ph, qqs, 5, 5, 5)


def test_rejects_nonuniform_source_coordinate(spherical_grid):
    """現行アルゴリズムが対応しない非一様球面座標を拒否する。"""
    rr, th, ph = spherical_grid
    rr = rr.copy()
    rr[2] += 0.1
    qqs = np.ones((rr.size, th.size, ph.size), dtype=np.float32)

    with pytest.raises(ValueError, match="uniformly spaced"):
        spherical2cartesian(rr, th, ph, qqs, 5, 5, 5)


def test_rejects_incomplete_phi_period(spherical_grid):
    """φ座標が完全な2π周期を覆わない入力を拒否する。"""
    rr, th, _ = spherical_grid
    ph = np.linspace(-0.5 * np.pi, 0.5 * np.pi, 9)
    qqs = np.ones((rr.size, th.size, ph.size), dtype=np.float32)

    with pytest.raises(ValueError, match=r"2\*pi"):
        spherical2cartesian(rr, th, ph, qqs, 5, 5, 5)


@pytest.mark.parametrize("sizes", [(1, 5, 5), (5, 1, 5), (5, 5, 1)])
def test_rejects_cartesian_axis_shorter_than_two(spherical_grid, sizes):
    """直交座標生成でゼロ除算になる2点未満の出力軸を拒否する。"""
    rr, th, ph = spherical_grid
    qqs = np.ones((rr.size, th.size, ph.size), dtype=np.float32)

    with pytest.raises(ValueError, match="size must be at least 2"):
        spherical2cartesian(rr, th, ph, qqs, *sizes)
