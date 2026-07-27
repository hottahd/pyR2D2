#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include "numpy_util.hpp"

namespace py = pybind11;

namespace interpolation_detail
{
using numpy_detail::CoordinateArray;

/*
 * 補間先の1座標について、補間元の下側格子点と線形補間重みを保持する。
 *
 * value = lower_weight * q[lower]
 *       + upper_weight * q[lower + 1]
 *
 * inside=falseは補間元領域の端点または領域外を表し、従来のFortran実装と
 * 同じく対応する出力をゼロにする。
 */
struct Bracket
{
  std::size_t lower = 0;
  double lower_weight = 0.0;
  double upper_weight = 0.0;
  bool inside = false;
};

/*
 * 補間元座標を検査し、各補間先座標の区間位置と重みを事前計算する。
 *
 * 現行Fortran版が前提としている、有限で狭義単調増加な補間元座標だけを
 * 受け付ける。補間先は並んでいる必要はないが、NaN・無限大は拒否する。
 * 補間元の両端と領域外はinside=falseにして、後段で出力をゼロにする。
 */
inline std::vector<Bracket> make_brackets(
    const CoordinateArray &source,
    const CoordinateArray &target,
    const char *source_name,
    const char *target_name)
{
  const py::ssize_t source_size = source.shape(0);
  const py::ssize_t target_size = target.shape(0);
  if (source_size < 2)
    throw py::value_error(
        std::string(source_name) + " must contain at least 2 points");

  const double *source_data = source.data();
  const double *target_data = target.data();
  if (!std::isfinite(source_data[0]))
    throw py::value_error(
        std::string(source_name) + " must be finite and strictly increasing");

  for (py::ssize_t i = 1; i < source_size; ++i)
  {
    if (!std::isfinite(source_data[i]) ||
        source_data[i] <= source_data[i - 1])
      throw py::value_error(
          std::string(source_name) + " must be finite and strictly increasing");
  }

  std::vector<Bracket> brackets(static_cast<std::size_t>(target_size));
  for (py::ssize_t i = 0; i < target_size; ++i)
  {
    const double value = target_data[i];
    if (!std::isfinite(value))
      throw py::value_error(std::string(target_name) + " must be finite");

    if (value <= source_data[0] || value >= source_data[source_size - 1])
      continue;

    const double *upper = std::upper_bound(
        source_data,
        source_data + source_size,
        value);
    const std::size_t lower =
        static_cast<std::size_t>(upper - source_data - 1);
    const double span = source_data[lower + 1] - source_data[lower];
    const double upper_weight = (value - source_data[lower]) / span;

    brackets[static_cast<std::size_t>(i)] = {
        lower,
        1.0 - upper_weight,
        upper_weight,
        true};
  }
  return brackets;
}

/*
 * C連続またはF連続な3次元配列を、直交する三つの1次元座標に沿って補間する本体。
 *
 * 入力qqはコピーせずに参照する。出力は補間先shape=(nxu, nyu, nzu)を持つ
 * C連続配列として新しく1個確保し、dtypeは入力と同じTにする。補間重みと
 * 8点の積和はdoubleで計算し、最後にTへ戻す。
 */
template <typename T>
py::array_t<T> interpolate_3d(
    const CoordinateArray &x,
    const CoordinateArray &y,
    const CoordinateArray &z,
    const CoordinateArray &xu,
    const CoordinateArray &yu,
    const CoordinateArray &zu,
    const py::array &qq)
{
  auto input = numpy_detail::borrow_contiguous_3d<T>(qq, "qq");

  const py::ssize_t nx = input.shape(0);
  const py::ssize_t ny = input.shape(1);
  const py::ssize_t nz = input.shape(2);
  if (x.shape(0) != nx || y.shape(0) != ny || z.shape(0) != nz)
    throw py::value_error(
        "source coordinate lengths must match the corresponding qq axes");

  const auto x_brackets = make_brackets(x, xu, "x", "xu");
  const auto y_brackets = make_brackets(y, yu, "y", "yu");
  const auto z_brackets = make_brackets(z, zu, "z", "zu");

  const py::ssize_t nxu = xu.shape(0);
  const py::ssize_t nyu = yu.shape(0);
  const py::ssize_t nzu = zu.shape(0);

  py::array_t<T> output({nxu, nyu, nzu});

  const T *input_data = input.data();
  T *output_data = output.mutable_data();
  const std::size_t input_stride_x =
      numpy_detail::element_stride(input, 0);
  const std::size_t input_stride_y =
      numpy_detail::element_stride(input, 1);
  const std::size_t input_stride_z =
      numpy_detail::element_stride(input, 2);
  const std::size_t output_stride_x =
      static_cast<std::size_t>(nyu) * static_cast<std::size_t>(nzu);
  const std::size_t output_stride_y = static_cast<std::size_t>(nzu);

  {
    /*
     * 各(i,j)が異なるz方向の出力列を担当するため書き込み競合はない。
     * 補間中はPython APIに触れないのでGILを解放する。
     */
    py::gil_scoped_release release;

#pragma omp parallel for collapse(2) schedule(static)
    for (py::ssize_t i = 0; i < nxu; ++i)
    {
      for (py::ssize_t j = 0; j < nyu; ++j)
      {
        const Bracket &bx = x_brackets[static_cast<std::size_t>(i)];
        const Bracket &by = y_brackets[static_cast<std::size_t>(j)];
        const std::size_t output_base =
            static_cast<std::size_t>(i) * output_stride_x +
            static_cast<std::size_t>(j) * output_stride_y;

        if (!bx.inside || !by.inside)
        {
          for (py::ssize_t k = 0; k < nzu; ++k)
            output_data[output_base + static_cast<std::size_t>(k)] = T(0);
          continue;
        }

        const std::size_t x0 = bx.lower;
        const std::size_t x1 = x0 + 1;
        const std::size_t y0 = by.lower;
        const std::size_t y1 = y0 + 1;

        for (py::ssize_t k = 0; k < nzu; ++k)
        {
          const Bracket &bz = z_brackets[static_cast<std::size_t>(k)];
          const std::size_t output_index =
              output_base + static_cast<std::size_t>(k);
          if (!bz.inside)
          {
            output_data[output_index] = T(0);
            continue;
          }

          const std::size_t z0 = bz.lower;
          const std::size_t z1 = z0 + 1;
          const auto index = [input_stride_x, input_stride_y, input_stride_z](
                                 const std::size_t ii,
                                 const std::size_t jj,
                                 const std::size_t kk)
          {
            return ii * input_stride_x +
                   jj * input_stride_y +
                   kk * input_stride_z;
          };

          const double value =
              static_cast<double>(input_data[index(x0, y0, z0)]) *
                  bx.lower_weight * by.lower_weight * bz.lower_weight +
              static_cast<double>(input_data[index(x0, y0, z1)]) *
                  bx.lower_weight * by.lower_weight * bz.upper_weight +
              static_cast<double>(input_data[index(x0, y1, z0)]) *
                  bx.lower_weight * by.upper_weight * bz.lower_weight +
              static_cast<double>(input_data[index(x0, y1, z1)]) *
                  bx.lower_weight * by.upper_weight * bz.upper_weight +
              static_cast<double>(input_data[index(x1, y0, z0)]) *
                  bx.upper_weight * by.lower_weight * bz.lower_weight +
              static_cast<double>(input_data[index(x1, y0, z1)]) *
                  bx.upper_weight * by.lower_weight * bz.upper_weight +
              static_cast<double>(input_data[index(x1, y1, z0)]) *
                  bx.upper_weight * by.upper_weight * bz.lower_weight +
              static_cast<double>(input_data[index(x1, y1, z1)]) *
                  bx.upper_weight * by.upper_weight * bz.upper_weight;
          output_data[output_index] = static_cast<T>(value);
        }
      }
    }
  }

  return output;
}

/*
 * Python引数を座標配列へ変換し、qqのdtypeに対応するテンプレート実体を選ぶ。
 *
 * 1次元座標は小さいためfloat64・C連続への変換を許すが、巨大なqqには
 * forcecastを使わない。対応外dtypeや逆エンディアンは暗黙コピーせず拒否する。
 */
inline py::array interpolation_dispatch(
    const py::array &x,
    const py::array &y,
    const py::array &z,
    const py::array &xu,
    const py::array &yu,
    const py::array &zu,
    const py::array &qq)
{
  const auto x_double = numpy_detail::ensure_coordinate_1d(x, "x");
  const auto y_double = numpy_detail::ensure_coordinate_1d(y, "y");
  const auto z_double = numpy_detail::ensure_coordinate_1d(z, "z");
  const auto xu_double = numpy_detail::ensure_coordinate_1d(xu, "xu");
  const auto yu_double = numpy_detail::ensure_coordinate_1d(yu, "yu");
  const auto zu_double = numpy_detail::ensure_coordinate_1d(zu, "zu");

  if (numpy_detail::is_native_float_dtype<float>(qq.dtype()))
    return interpolate_3d<float>(
        x_double,
        y_double,
        z_double,
        xu_double,
        yu_double,
        zu_double,
        qq);
  if (numpy_detail::is_native_float_dtype<double>(qq.dtype()))
    return interpolate_3d<double>(
        x_double,
        y_double,
        z_double,
        xu_double,
        yu_double,
        zu_double,
        qq);

  throw py::type_error(
      "qq must have native-endian dtype float32 or float64 (" +
      numpy_detail::dtype_description(qq.dtype()) + ")");
}
} // namespace interpolation_detail
