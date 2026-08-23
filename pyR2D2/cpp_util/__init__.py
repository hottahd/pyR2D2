"""Utilities backed by the C++ extension module."""

from .cpp_util import (
    EOS,
    YinYang,
    d_x as _d_x_3d,
    d_y as _d_y_3d,
    d_z as _d_z_3d,
    eval_tau,
    interp,
    spherical2cartesian,
    trace_field_line,
    trace_field_lines,
    vertical_upward_rte,
)


def _derivative_input_3d(coordinate, qq, axis):
    """1D・2D入力へsingleton軸を加え、C++コア用の3D viewを返す。"""
    if qq.ndim == 3:
        return qq

    if qq.ndim == 1:
        if axis == 0:
            return qq[:, None, None]
        if axis == 1:
            return qq[None, :, None]
        return qq[None, None, :]

    if qq.ndim == 2:
        if axis == 0:
            return qq[:, :, None]
        if axis == 2:
            return qq[:, None, :]

        # d_yでは(y, z)と(x, y)の二通りがある。正方配列など両方に
        # 一致する場合は、従来どおり第0軸をyと解釈する。
        coordinate_size = len(coordinate)
        if qq.shape[0] == coordinate_size:
            return qq[None, :, :]
        if qq.shape[1] == coordinate_size:
            return qq[:, :, None]
        raise ValueError(
            "coordinate length must match one axis of a 2D input for d_y"
        )

    raise ValueError("qq must be a 1D, 2D, or 3D array")


def _differentiate(coordinate, qq, axis, derivative_3d):
    """3D C++コアを呼び、追加コピーなしで入力shapeのviewへ戻す。"""
    input_3d = _derivative_input_3d(coordinate, qq, axis)
    result_3d = derivative_3d(coordinate, input_3d)
    if result_3d.shape == qq.shape:
        return result_3d
    return result_3d.reshape(qq.shape)


def d_x(x, qq):
    """Differentiate an array along its x axis.

    Parameters
    ----------
    x : numpy.ndarray
        One-dimensional x-coordinate array.
    qq : numpy.ndarray
        One-, two-, or three-dimensional, C- or F-contiguous,
        native-endian float32 or float64 array. Its first axis is x.

    Returns
    -------
    numpy.ndarray
        C-contiguous derivative with the same shape and dtype as ``qq``.

    Notes
    -----
    Singleton axes are added only as views while the three-dimensional C++
    core runs. No input copy is made. The first and last two x cells are zero.
    """
    return _differentiate(x, qq, axis=0, derivative_3d=_d_x_3d)


def d_y(y, qq):
    """Differentiate an array along its y axis.

    Parameters
    ----------
    y : numpy.ndarray
        One-dimensional y-coordinate array.
    qq : numpy.ndarray
        One-, two-, or three-dimensional, C- or F-contiguous,
        native-endian float32 or float64 array. A 2D input is interpreted as
        ``(y, z)`` if its first axis matches ``y``; otherwise as ``(x, y)``
        if its second axis matches ``y``.

    Returns
    -------
    numpy.ndarray
        C-contiguous derivative with the same shape and dtype as ``qq``.

    Notes
    -----
    Singleton axes are added only as views while the three-dimensional C++
    core runs. No input copy is made. The first and last two y cells are zero.
    """
    return _differentiate(y, qq, axis=1, derivative_3d=_d_y_3d)


def d_z(z, qq):
    """Differentiate an array along its z axis.

    Parameters
    ----------
    z : numpy.ndarray
        One-dimensional z-coordinate array.
    qq : numpy.ndarray
        One-, two-, or three-dimensional, C- or F-contiguous,
        native-endian float32 or float64 array. Its last axis is z.

    Returns
    -------
    numpy.ndarray
        C-contiguous derivative with the same shape and dtype as ``qq``.

    Notes
    -----
    Singleton axes are added only as views while the three-dimensional C++
    core runs. No input copy is made. The first and last two z cells are zero.
    """
    return _differentiate(z, qq, axis=2, derivative_3d=_d_z_3d)


__all__ = [
    "EOS",
    "YinYang",
    "d_x",
    "d_y",
    "d_z",
    "eval_tau",
    "interp",
    "spherical2cartesian",
    "vertical_upward_rte",
    "trace_field_line",
    "trace_field_lines",
]
