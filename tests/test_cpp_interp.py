"""C++三線形補間器の精度、dtype、および入力条件を確認するテスト。"""

import numpy as np
import pytest

from pyR2D2.cpp_util import interp


@pytest.fixture
def grids():
    """非一様な補間元座標と、その内部にある補間先座標を返す。"""
    x = np.array([0.0, 0.7, 2.0, 4.0])
    y = np.array([-2.0, -0.5, 1.5])
    z = np.array([1.0, 1.8, 3.5, 6.0])
    xu = np.array([0.2, 1.0, 3.0])
    yu = np.array([-1.5, 0.5, 1.0])
    zu = np.array([1.2, 2.5, 5.0])
    return (x, y, z), (xu, yu, zu)


def trilinear_function(x, y, z):
    """三線形補間で厳密に再現できる試験関数を返す。"""
    return (
        1.0
        + 2.0 * x[:, None, None]
        - 3.0 * y[None, :, None]
        + 0.5 * z[None, None, :]
        + 0.25 * x[:, None, None] * y[None, :, None] * z[None, None, :]
    )


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
def test_nonuniform_trilinear_interpolation_preserves_dtype(grids, dtype):
    """非一様格子を補間し、入力と同じdtypeのC連続配列を返すことを確認する。"""
    (x, y, z), (xu, yu, zu) = grids
    qq = np.ascontiguousarray(trilinear_function(x, y, z), dtype=dtype)

    result = interp(x, y, z, xu, yu, zu, qq)
    expected = trilinear_function(xu, yu, zu)

    assert result.dtype == dtype
    assert result.flags.c_contiguous
    tolerance = 2.0e-6 if dtype == np.float32 else 1.0e-12
    np.testing.assert_allclose(result, expected, rtol=0.0, atol=tolerance)


@pytest.mark.parametrize("axis", [0, 1, 2])
def test_source_boundaries_and_outside_are_zero(grids, axis):
    """補間元領域の端点および領域外に対応する出力面がゼロになることを確認する。"""
    (x, y, z), targets = grids
    coordinates = (x, y, z)
    qq = np.ascontiguousarray(trilinear_function(x, y, z), dtype=np.float32)

    target_coordinates = list(targets)
    source = coordinates[axis]
    target_coordinates[axis] = np.array(
        [
            source[0] - 0.1,
            source[0],
            0.5 * (source[0] + source[1]),
            source[-1],
            source[-1] + 0.1,
        ]
    )

    result = interp(x, y, z, *target_coordinates, qq)

    for index in (0, 1, 3, 4):
        np.testing.assert_array_equal(np.take(result, index, axis=axis), 0.0)

    assert np.any(np.take(result, 2, axis=axis) != 0.0)


@pytest.mark.parametrize("invalid_dtype", [np.int32, np.int64, np.float16])
def test_rejects_unsupported_input_dtype(grids, invalid_dtype):
    """暗黙の巨大変換を避けるため、float32/float64以外を拒否する。"""
    (x, y, z), (xu, yu, zu) = grids
    qq = np.ones((x.size, y.size, z.size), dtype=invalid_dtype)

    with pytest.raises(TypeError, match="float32 or float64"):
        interp(x, y, z, xu, yu, zu, qq)


def test_accepts_f_contiguous_input(grids):
    """F連続入力をコピーせず読み、C連続入力と同じ結果を返す。"""
    (x, y, z), (xu, yu, zu) = grids
    qq_c = np.ascontiguousarray(trilinear_function(x, y, z), dtype=np.float32)
    qq_f = np.asfortranarray(qq_c)

    expected = interp(x, y, z, xu, yu, zu, qq_c)
    result = interp(x, y, z, xu, yu, zu, qq_f)

    assert qq_f.flags.f_contiguous
    np.testing.assert_array_equal(result, expected)


def test_rejects_noncontiguous_input(grids):
    """C連続でもF連続でもない飛び飛びviewを拒否する。"""
    (x, y, z), (xu, yu, zu) = grids
    storage = np.ones((x.size, y.size, 2 * z.size), dtype=np.float32)
    qq = storage[:, :, ::2]

    assert not qq.flags.c_contiguous
    assert not qq.flags.f_contiguous
    with pytest.raises(ValueError, match="C- or F-contiguous"):
        interp(x, y, z, xu, yu, zu, qq)


def test_accepts_native_float_dtype_with_metadata(grids):
    """dtype metadataを持つR2D2配列もコピーせず受け付ける。"""
    (x, y, z), (xu, yu, zu) = grids
    dtype = np.dtype(np.float32, metadata={"source": "r2d2"})
    qq = np.ones((x.size, y.size, z.size), dtype=dtype)

    result = interp(x, y, z, xu, yu, zu, qq)

    assert result.dtype == np.float32


def test_rejects_non_3d_input(grids):
    """補間対象は3次元配列だけを受け付ける。"""
    (x, y, z), (xu, yu, zu) = grids
    qq = np.ones((x.size, y.size), dtype=np.float32)

    with pytest.raises(ValueError, match="3D"):
        interp(x, y, z, xu, yu, zu, qq)


def test_rejects_source_coordinate_length_mismatch(grids):
    """補間元座標長とqqの対応軸長が異なる入力を拒否する。"""
    (x, y, z), (xu, yu, zu) = grids
    qq = np.ones((x.size, y.size, z.size), dtype=np.float32)

    with pytest.raises(ValueError, match="source coordinate lengths"):
        interp(x[:-1], y, z, xu, yu, zu, qq)


@pytest.mark.parametrize(
    "invalid_x",
    [
        np.array([0.0]),
        np.array([0.0, 1.0, 1.0, 2.0]),
        np.array([0.0, 2.0, 1.0, 3.0]),
        np.array([0.0, 1.0, np.nan, 3.0]),
    ],
)
def test_rejects_invalid_source_coordinate(grids, invalid_x):
    """短すぎる、重複、非単調、非有限な補間元座標を拒否する。"""
    (_, y, z), (xu, yu, zu) = grids
    qq = np.ones((invalid_x.size, y.size, z.size), dtype=np.float32)

    with pytest.raises(ValueError, match="at least 2|finite and strictly increasing"):
        interp(invalid_x, y, z, xu, yu, zu, qq)


def test_rejects_nonfinite_target_coordinate(grids):
    """補間先座標のNaN・無限大を拒否する。"""
    (x, y, z), (_, yu, zu) = grids
    xu = np.array([0.2, np.nan, 3.0])
    qq = np.ones((x.size, y.size, z.size), dtype=np.float32)

    with pytest.raises(ValueError, match="xu must be finite"):
        interp(x, y, z, xu, yu, zu, qq)
