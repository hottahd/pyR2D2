"""C++微分器の数値精度、dtype、および入力条件を確認するテスト。"""

import numpy as np
import pytest

from pyR2D2.cpp_util import d_x, d_y, d_z


@pytest.fixture
def nonuniform_grid():
    """各軸が5点以上ある、小さな非一様3次元格子を返す。"""
    x = np.array([-2.0, -1.4, -0.7, 0.1, 1.0, 2.2, 3.7])
    y = np.array([-1.8, -1.1, -0.3, 0.6, 1.6, 2.9, 4.5, 6.4])
    z = np.array([-2.5, -1.7, -0.8, 0.2, 1.3, 2.5, 3.8, 5.2, 6.7])
    return x, y, z


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
@pytest.mark.parametrize(
    ("derivative", "axis", "factor"),
    [
        (d_x, 0, 1.0),
        (d_y, 1, 2.0),
        (d_z, 2, -0.5),
    ],
)
def test_cubic_on_nonuniform_grid_preserves_dtype(
    nonuniform_grid,
    dtype,
    derivative,
    axis,
    factor,
):
    """非一様格子の3次式を微分し、入力と同じdtypeで返すことを確認する。"""
    x, y, z = nonuniform_grid
    coordinates = (x, y, z)
    q = np.ascontiguousarray(
        x[:, None, None] ** 3
        + 2.0 * y[None, :, None] ** 3
        - 0.5 * z[None, None, :] ** 3,
        dtype=dtype,
    )

    result = derivative(coordinates[axis], q)

    expected = np.zeros(q.shape, dtype=np.float64)
    derivative_shape = [1, 1, 1]
    derivative_shape[axis] = coordinates[axis].size
    analytic = (3.0 * factor * coordinates[axis] ** 2).reshape(derivative_shape)

    interior = [slice(None)] * 3
    interior[axis] = slice(2, -2)
    expected[tuple(interior)] = np.broadcast_to(analytic, q.shape)[tuple(interior)]

    assert result.dtype == dtype
    assert result.flags.c_contiguous
    tolerance = 3.0e-5 if dtype == np.float32 else 1.0e-12
    np.testing.assert_allclose(result, expected, rtol=0.0, atol=tolerance)


@pytest.mark.parametrize(
    ("derivative", "axis"),
    [
        (d_x, 0),
        (d_y, 1),
        (d_z, 2),
    ],
)
def test_two_boundary_cells_are_zero(nonuniform_grid, derivative, axis):
    """微分軸の両端2格子がゼロになることを確認する。"""
    x, y, z = nonuniform_grid
    coordinates = (x, y, z)
    q = np.ascontiguousarray(
        x[:, None, None] + y[None, :, None] + z[None, None, :],
        dtype=np.float32,
    )

    result = derivative(coordinates[axis], q)

    boundary = [slice(None)] * 3
    for index in (0, 1, -2, -1):
        boundary[axis] = index
        np.testing.assert_array_equal(result[tuple(boundary)], 0.0)


@pytest.mark.parametrize("invalid_dtype", [np.int32, np.int64, np.float16])
def test_rejects_unsupported_input_dtype(nonuniform_grid, invalid_dtype):
    """暗黙の巨大変換を避けるため、float32/float64以外を拒否する。"""
    x, y, z = nonuniform_grid
    q = np.ones((x.size, y.size, z.size), dtype=invalid_dtype)

    with pytest.raises(TypeError, match="float32 or float64"):
        d_x(x, q)


def test_rejects_non_c_contiguous_input(nonuniform_grid):
    """暗黙コピーを避けるため、C連続でない入力を拒否する。"""
    x, y, z = nonuniform_grid
    q = np.ones((x.size, y.size, z.size), dtype=np.float32, order="F")

    with pytest.raises(ValueError, match="C-contiguous"):
        d_x(x, q)


def test_accepts_native_float_dtype_with_metadata(nonuniform_grid):
    """dtypeの付随情報だけを理由に巨大な入力配列を拒否しない。"""
    x, y, z = nonuniform_grid
    dtype = np.dtype(np.float32, metadata={"source": "r2d2"})
    q = np.ones((x.size, y.size, z.size), dtype=dtype)

    result = d_x(x, q)

    assert result.dtype == np.float32


def test_rejects_non_3d_input(nonuniform_grid):
    """今回の初期実装は3次元配列だけを対象とする。"""
    x, _, _ = nonuniform_grid
    q = np.ones((x.size, 3), dtype=np.float32)

    with pytest.raises(ValueError, match="3D"):
        d_x(x, q)


def test_rejects_coordinate_length_mismatch(nonuniform_grid):
    """座標長と微分対象軸の長さが異なる入力を拒否する。"""
    x, y, z = nonuniform_grid
    q = np.ones((x.size, y.size, z.size), dtype=np.float32)

    with pytest.raises(ValueError, match="coordinate length"):
        d_x(x[:-1], q)


def test_rejects_short_differentiated_axis():
    """4点公式を適用できない5点未満の軸を拒否する。"""
    x = np.linspace(0.0, 1.0, 4)
    q = np.ones((x.size, 5, 6), dtype=np.float32)

    with pytest.raises(ValueError, match="at least 5"):
        d_x(x, q)


@pytest.mark.parametrize(
    "coordinate",
    [
        np.array([0.0, 1.0, 1.0, 2.0, 3.0]),
        np.array([0.0, 1.0, 0.5, 2.0, 3.0]),
        np.array([0.0, 1.0, np.nan, 2.0, 3.0]),
    ],
)
def test_rejects_invalid_coordinate(coordinate):
    """重複、非単調、非有限な座標を拒否する。"""
    q = np.ones((coordinate.size, 5, 6), dtype=np.float32)

    with pytest.raises(ValueError, match="finite and strictly monotonic"):
        d_x(coordinate, q)
