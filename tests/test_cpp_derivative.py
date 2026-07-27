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


@pytest.mark.parametrize(
    ("derivative", "axis"),
    [(d_x, 0), (d_y, 1), (d_z, 2)],
)
def test_accepts_f_contiguous_input(nonuniform_grid, derivative, axis):
    """F連続入力をコピーせず読み、C連続入力と同じ結果を返す。"""
    x, y, z = nonuniform_grid
    coordinates = (x, y, z)
    q_c = np.ascontiguousarray(
        x[:, None, None] ** 3
        + 2.0 * y[None, :, None] ** 3
        - 0.5 * z[None, None, :] ** 3,
        dtype=np.float32,
    )
    q_f = np.asfortranarray(q_c)

    expected = derivative(coordinates[axis], q_c)
    result = derivative(coordinates[axis], q_f)

    assert q_f.flags.f_contiguous
    np.testing.assert_array_equal(result, expected)


def test_rejects_noncontiguous_input(nonuniform_grid):
    """C連続でもF連続でもない飛び飛びviewを拒否する。"""
    x, y, z = nonuniform_grid
    storage = np.ones((x.size, y.size, 2 * z.size), dtype=np.float32)
    q = storage[:, :, ::2]

    assert not q.flags.c_contiguous
    assert not q.flags.f_contiguous
    with pytest.raises(ValueError, match="C- or F-contiguous"):
        d_x(x, q)


def test_accepts_native_float_dtype_with_metadata(nonuniform_grid):
    """dtypeの付随情報だけを理由に巨大な入力配列を拒否しない。"""
    x, y, z = nonuniform_grid
    dtype = np.dtype(np.float32, metadata={"source": "r2d2"})
    q = np.ones((x.size, y.size, z.size), dtype=dtype)

    result = d_x(x, q)

    assert result.dtype == np.float32


@pytest.mark.parametrize(
    ("derivative", "coordinate_index"),
    [(d_x, 0), (d_y, 1), (d_z, 2)],
)
def test_1d_input_returns_1d(nonuniform_grid, derivative, coordinate_index):
    """1D入力を同じshape・dtypeの1D配列として返す。"""
    coordinate = nonuniform_grid[coordinate_index]
    q = np.asarray(coordinate**3, dtype=np.float32)

    result = derivative(coordinate, q)

    expected = np.zeros_like(q)
    expected[2:-2] = 3.0 * coordinate[2:-2] ** 2
    assert result.shape == q.shape
    assert result.ndim == 1
    assert result.dtype == q.dtype
    np.testing.assert_allclose(result, expected, rtol=0.0, atol=3.0e-5)


@pytest.mark.parametrize(
    ("derivative", "coordinate_index", "input_shape"),
    [
        (d_x, 0, (7, 5)),
        (d_y, 1, (8, 5)),  # (y, z)
        (d_y, 1, (6, 8)),  # (x, y)
        (d_z, 2, (5, 9)),
    ],
)
def test_2d_input_returns_2d(
    nonuniform_grid, derivative, coordinate_index, input_shape
):
    """対応する微分軸を判定し、2D入力と同じshapeで返す。"""
    coordinate = nonuniform_grid[coordinate_index]
    q = np.asfortranarray(
        np.arange(np.prod(input_shape), dtype=np.float64).reshape(input_shape)
    )

    result = derivative(coordinate, q)

    assert result.shape == q.shape
    assert result.ndim == 2
    assert result.dtype == q.dtype
    assert result.flags.c_contiguous


def test_d_y_square_2d_input_prefers_first_axis():
    """両軸長がyと一致する場合は第0軸をyとして微分する。"""
    y = np.linspace(-1.0, 1.0, 6)
    q = np.broadcast_to(y[:, None] ** 3, (6, 6)).copy()

    result = d_y(y, q)

    expected = np.zeros_like(q)
    expected[2:-2, :] = 3.0 * y[2:-2, None] ** 2
    np.testing.assert_allclose(result, expected, rtol=0.0, atol=1.0e-12)


def test_rejects_dimensions_outside_supported_range(nonuniform_grid):
    """0Dおよび4D配列を明示的に拒否する。"""
    x, _, _ = nonuniform_grid
    for q in (np.array(1.0), np.ones((x.size, 1, 1, 1))):
        with pytest.raises(ValueError, match="1D, 2D, or 3D"):
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
