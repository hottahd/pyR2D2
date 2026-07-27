"""現行Fortran微分の挙動を、C++移植前の基準として固定するテスト。"""

import numpy as np
import pytest

from pyR2D2.fortran_util import d_x, d_y, d_z


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
def test_cubic_on_nonuniform_grid_matches_analytic_derivative(
    nonuniform_grid,
    dtype,
    derivative,
    axis,
    factor,
):
    """4点公式が非一様格子上の3次多項式を微分できることを確認する。

    現行Fortranでは微分係数rk1--rk4がreal(4)なので、float64入力でも
    丸め誤差が残る。ここでは移植前の挙動を固定するため、その誤差を
    含めて十分に厳しい絶対許容値を設定する。
    """
    x, y, z = nonuniform_grid
    coordinates = (x, y, z)

    q = (
        x[:, None, None] ** 3
        + 2.0 * y[None, :, None] ** 3
        - 0.5 * z[None, None, :] ** 3
    ).astype(dtype)

    result = derivative(coordinates[axis], q)

    expected = np.zeros_like(result)
    derivative_shape = [1, 1, 1]
    derivative_shape[axis] = coordinates[axis].size
    analytic = (3.0 * factor * coordinates[axis] ** 2).reshape(derivative_shape)

    interior = [slice(None)] * 3
    interior[axis] = slice(2, -2)
    expected[tuple(interior)] = np.broadcast_to(analytic, result.shape)[
        tuple(interior)
    ]

    assert result.dtype == np.float64
    np.testing.assert_allclose(result, expected, rtol=0.0, atol=3.0e-5)


@pytest.mark.parametrize(
    ("derivative", "axis"),
    [
        (d_x, 0),
        (d_y, 1),
        (d_z, 2),
    ],
)
def test_two_boundary_cells_are_zero(nonuniform_grid, derivative, axis):
    """現行実装が微分軸の両端2格子をゼロのまま返すことを確認する。"""
    x, y, z = nonuniform_grid
    coordinates = (x, y, z)
    q = (
        x[:, None, None]
        + y[None, :, None]
        + z[None, None, :]
    ).astype(np.float32)

    result = derivative(coordinates[axis], q)

    boundary = [slice(None)] * 3
    for index in (0, 1, -2, -1):
        boundary[axis] = index
        np.testing.assert_array_equal(result[tuple(boundary)], 0.0)
