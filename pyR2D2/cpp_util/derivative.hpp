#pragma once

#include <cmath>
#include <cstddef>
#include <vector>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include "numpy_util.hpp"

namespace py = pybind11;

namespace derivative_detail
{
using numpy_detail::CoordinateArray;

constexpr py::ssize_t margin = 2;

struct Weights
{
  double rk1;
  double rk2;
  double rk3;
  double rk4;
};

/*
 * 1次元の非一様座標から、各内部格子点の4点中心微分係数を事前計算する。
 *
 * 微分点iそのものの値は使わず、i+2、i+1、i-1、i-2の値に掛ける4係数を
 * Weightsとして返す。係数は座標精度と積和精度を保つため常にdoubleで計算する。
 *
 * coordinateは小さな1次元配列なので、呼び出し側で必要に応じてC連続の
 * float64配列へ変換することを許している。一方、重複座標、非単調座標、
 * NaN・無限大、および公式を適用できない5点未満の座標はここで拒否する。
 *
 * 両端margin点の係数は実際の微分では使わず、出力値をゼロにする。
 */
inline std::vector<Weights> make_weights(
    const CoordinateArray &coordinate)
{
  if (coordinate.ndim() != 1)
    throw py::value_error("coordinate must be a 1D array");

  const py::ssize_t size = coordinate.shape(0);
  if (size < 2 * margin + 1)
    throw py::value_error("the differentiated axis must contain at least 5 points");

  const double *x = coordinate.data();
  const double first_delta = x[1] - x[0];
  if (!std::isfinite(first_delta) || first_delta == 0.0)
    throw py::value_error("coordinate must be finite and strictly monotonic");

  const bool increasing = first_delta > 0.0;
  for (py::ssize_t i = 1; i < size; ++i)
  {
    const double delta = x[i] - x[i - 1];
    if (!std::isfinite(x[i]) || delta == 0.0 || (delta > 0.0) != increasing)
      throw py::value_error("coordinate must be finite and strictly monotonic");
  }

  std::vector<Weights> weights(static_cast<std::size_t>(size));
  for (py::ssize_t i = margin; i < size - margin; ++i)
  {
    const double a = x[i + 2] - x[i];
    const double b = x[i + 1] - x[i];
    const double c = x[i] - x[i - 1];
    const double d = x[i] - x[i - 2];

    weights[static_cast<std::size_t>(i)] = {
        (c * d - b * (c + d)) / (a - b) / (a + c) / (a + d),
        (a * (c + d) - c * d) / (a - b) / (b + c) / (b + d),
        (b * d + a * (d - b)) / (a + c) / (b + c) / (c - d),
        (a * (b - c) - b * c) / (c - d) / (a + d) / (b + d)};
  }
  return weights;
}

/*
 * 3次元配列を1次元化したバッファ上で、中央の1点だけを微分する。
 *
 * centerは微分対象点の一次元添字、strideは微分方向へ1格子進むための
 * 要素数である。C連続・F連続の違いは呼び出し側がNumPyの実strideを
 * 渡すことで吸収し、この関数はどちらのメモリ配置にも同じ式を使う。
 *
 * 入力がfloat32でも係数との積和はdoubleで行い、最後にTへ戻す。
 * したがって出力dtypeは入力dtypeと同じになる。
 */
template <typename T>
inline T evaluate(
    const T *input,
    const std::size_t center,
    const std::size_t stride,
    const Weights &weight)
{
  const double value =
      weight.rk1 * static_cast<double>(input[center + 2 * stride]) +
      weight.rk2 * static_cast<double>(input[center + stride]) +
      weight.rk3 * static_cast<double>(input[center - stride]) +
      weight.rk4 * static_cast<double>(input[center - 2 * stride]);
  return static_cast<T>(value);
}

/*
 * C連続またはF連続な3次元NumPy配列を、指定された1軸方向へ微分する本体。
 *
 * Tは入出力の要素型（floatまたはdouble）、Axisは微分軸
 * （0=x、1=y、2=z）である。Axisはコンパイル時定数なので、if constexprで
 * 選ばれなかった軸のコードは生成されない。
 *
 * 入力qqはコピーせずに参照する。出力についてはqqと同じshape・dtypeの
 * C連続配列を新しく1個確保する。微分軸の両端margin=2格子は、従来の
 * Fortran実装と同じくゼロにする。
 */
template <typename T, std::size_t Axis>
py::array_t<T> derivative_3d(
    const CoordinateArray &coordinate,
    const py::array &qq)
{
  static_assert(Axis < 3);

  auto input = numpy_detail::borrow_contiguous_3d<T>(qq, "qq");

  const py::ssize_t nx = input.shape(0);
  const py::ssize_t ny = input.shape(1);
  const py::ssize_t nz = input.shape(2);
  const py::ssize_t axis_size = Axis == 0 ? nx : (Axis == 1 ? ny : nz);

  if (coordinate.shape(0) != axis_size)
    throw py::value_error("coordinate length must match the differentiated axis");

  const auto weights = make_weights(coordinate);

  // 入力と同じshape・dtypeの出力を新規確保する。入力データは複製しない。
  py::array_t<T> output({nx, ny, nz});

  const T *input_data = input.data();
  T *output_data = output.mutable_data();

  /*
   * 入力はC連続またはF連続なので、NumPyが持つ実strideを要素数単位で使う。
   * 出力は常にC連続であり、入力と出力の一次元添字を分けて計算する。
   */
  const std::size_t input_stride_i =
      numpy_detail::element_stride(input, 0);
  const std::size_t input_stride_j =
      numpy_detail::element_stride(input, 1);
  const std::size_t input_stride_k =
      numpy_detail::element_stride(input, 2);
  const std::size_t output_stride_i =
      static_cast<std::size_t>(ny) * static_cast<std::size_t>(nz);
  const std::size_t output_stride_j = static_cast<std::size_t>(nz);

  {
    /*
     * 以下ではPython APIに触れず、生ポインタだけを操作するためGILを解放する。
     * (i,j)の組をOpenMPスレッドへ静的に分配する。各組は異なるz方向の列を
     * 書き込むので、複数スレッドが同じ出力要素へ書く競合は発生しない。
     */
    py::gil_scoped_release release;

#pragma omp parallel for collapse(2) schedule(static)
    for (py::ssize_t i = 0; i < nx; ++i)
    {
      for (py::ssize_t j = 0; j < ny; ++j)
      {
        const std::size_t output_base =
            static_cast<std::size_t>(i) * output_stride_i +
            static_cast<std::size_t>(j) * output_stride_j;
        const std::size_t input_base =
            static_cast<std::size_t>(i) * input_stride_i +
            static_cast<std::size_t>(j) * input_stride_j;

        if constexpr (Axis == 0)
        {
          if (i < margin || i >= nx - margin)
          {
            for (py::ssize_t k = 0; k < nz; ++k)
              output_data[output_base + static_cast<std::size_t>(k)] = T(0);
            continue;
          }

          const Weights &weight = weights[static_cast<std::size_t>(i)];
          for (py::ssize_t k = 0; k < nz; ++k)
          {
            const std::size_t input_center =
                input_base + static_cast<std::size_t>(k) * input_stride_k;
            const std::size_t output_index =
                output_base + static_cast<std::size_t>(k);
            output_data[output_index] = evaluate(
                input_data, input_center, input_stride_i, weight);
          }
        }
        else if constexpr (Axis == 1)
        {
          if (j < margin || j >= ny - margin)
          {
            for (py::ssize_t k = 0; k < nz; ++k)
              output_data[output_base + static_cast<std::size_t>(k)] = T(0);
            continue;
          }

          const Weights &weight = weights[static_cast<std::size_t>(j)];
          for (py::ssize_t k = 0; k < nz; ++k)
          {
            const std::size_t input_center =
                input_base + static_cast<std::size_t>(k) * input_stride_k;
            const std::size_t output_index =
                output_base + static_cast<std::size_t>(k);
            output_data[output_index] = evaluate(
                input_data, input_center, input_stride_j, weight);
          }
        }
        else
        {
          output_data[output_base] = T(0);
          output_data[output_base + 1] = T(0);
          output_data[output_base + static_cast<std::size_t>(nz - 2)] = T(0);
          output_data[output_base + static_cast<std::size_t>(nz - 1)] = T(0);

          for (py::ssize_t k = margin; k < nz - margin; ++k)
          {
            const std::size_t input_center =
                input_base + static_cast<std::size_t>(k) * input_stride_k;
            const std::size_t output_index =
                output_base + static_cast<std::size_t>(k);
            output_data[output_index] = evaluate(
                input_data,
                input_center,
                input_stride_k,
                weights[static_cast<std::size_t>(k)]);
          }
        }
      }
    }
  }

  return output;
}

/*
 * Pythonから受け取ったNumPy dtypeを、対応するC++テンプレート実体へ振り分ける。
 *
 * coordinateは軸長程度の小配列なので、必要ならC連続float64へ変換する。
 * 巨大になり得るqqにはforcecastを使わず、ネイティブfloat32なら
 * derivative_3d<float, Axis>、float64ならderivative_3d<double, Axis>を呼ぶ。
 * 対応外dtypeは暗黙変換・暗黙コピーを行わず、診断情報を付けて拒否する。
 */
template <std::size_t Axis>
py::array derivative_dispatch(const py::array &coordinate, const py::array &qq)
{
  const auto coordinate_double =
      numpy_detail::ensure_coordinate_1d(coordinate, "coordinate");

  if (numpy_detail::is_native_float_dtype<float>(qq.dtype()))
    return derivative_3d<float, Axis>(coordinate_double, qq);
  if (numpy_detail::is_native_float_dtype<double>(qq.dtype()))
    return derivative_3d<double, Axis>(coordinate_double, qq);

  throw py::type_error(
      "qq must have native-endian dtype float32 or float64 (" +
      numpy_detail::dtype_description(qq.dtype()) + ")");
}
} // namespace derivative_detail
