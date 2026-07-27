#pragma once

#include <bit>
#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace derivative_detail
{
constexpr py::ssize_t margin = 2;

struct Weights
{
  double rk1;
  double rk2;
  double rk3;
  double rk4;
};

/*
 * NumPyのdtypeが、C++の浮動小数点型Tとしてコピーなしで安全に読めるか判定する。
 *
 * dtypeの完全一致（dtype.is(py::dtype::of<T>())）は使用しない。R2D2から
 * 読み込んだ配列のようにdtype metadataを持つ配列は、実体が通常のfloat32
 * またはfloat64でも、metadataを持たない標準dtypeとは別物と判定されるためである。
 * ここでは実際のメモリを読むために必要な次の三条件だけを確認する。
 *
 *   1. kindが `f`、すなわち浮動小数点型である。
 *   2. 1要素のバイト数がC++のT（floatまたはdouble）と一致する。
 *   3. 配列のバイト順が現在のCPUと一致する。
 *
 * NumPyのdtype.byteorder()が返す主な文字は次の通り。
 *
 *   `=` : CPUネイティブのバイト順という相対的な指定
 *   `<` : リトルエンディアンという絶対的な指定
 *   `>` : ビッグエンディアンという絶対的な指定
 *   `|` : バイト順を考える必要がない型
 *
 * 例えばリトルエンディアンCPUでは、`=`と`<`はいずれも実際には
 * リトルエンディアンを表す。ただしNumPyから返される文字表現は異なり得るため、
 * `=`を直接許可する条件と、`<`がCPUのstd::endian::littleと一致する条件の
 * 両方が必要になる。ビッグエンディアンについても同様に扱う。
 *
 * CPUと逆のバイト順を許可すると、入力をT*として参照したときに各要素の
 * バイト列を誤って解釈する。ここでは巨大な暗黙byteswapコピーを避けるため、
 * そのような配列は変換せず明示的に拒否する。
 */
template <typename T>
inline bool is_native_float_dtype(const py::dtype &dtype)
{
  const char byte_order = dtype.byteorder();
  const bool native_byte_order =
      byte_order == '=' ||
      byte_order == '|' ||
      (byte_order == '<' && std::endian::native == std::endian::little) ||
      (byte_order == '>' && std::endian::native == std::endian::big);
  return dtype.kind() == 'f' &&
         dtype.itemsize() == static_cast<py::ssize_t>(sizeof(T)) &&
         native_byte_order;
}

/*
 * dtype検査で例外を送出するときの補助関数。
 *
 * Python側から渡されたdtypeについて、型名だけでなくkind、1要素のバイト数、
 * byteorderも一つの文字列にまとめる。計算には使わず、対応外dtypeの原因を
 * エラーメッセージから判別しやすくするためだけに使用する。
 */
inline std::string dtype_description(const py::dtype &dtype)
{
  return "dtype=" + py::str(dtype).cast<std::string>() +
         ", kind=" + std::string(1, dtype.kind()) +
         ", itemsize=" + std::to_string(dtype.itemsize()) +
         ", byteorder=" + std::string(1, dtype.byteorder());
}

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
    const py::array_t<double, py::array::c_style | py::array::forcecast> &coordinate)
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
 * 要素数である。C連続なshape=(nx, ny, nz)の配列では、x、y、z方向の
 * strideはそれぞれny*nz、nz、1となる。
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
 * C連続な3次元NumPy配列全体を、指定された1軸方向へ微分する本体。
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
    const py::array_t<double, py::array::c_style | py::array::forcecast> &coordinate,
    const py::array &qq)
{
  static_assert(Axis < 3);

  if (qq.ndim() != 3)
    throw py::value_error("qq must be a 3D array");
  if (!(qq.flags() & py::array::c_style))
    throw py::value_error("qq must be C-contiguous");
  if (!is_native_float_dtype<T>(qq.dtype()))
    throw py::type_error("qq dtype does not match the selected derivative kernel");

  const py::ssize_t nx = qq.shape(0);
  const py::ssize_t ny = qq.shape(1);
  const py::ssize_t nz = qq.shape(2);
  const py::ssize_t axis_size = Axis == 0 ? nx : (Axis == 1 ? ny : nz);

  if (coordinate.shape(0) != axis_size)
    throw py::value_error("coordinate length must match the differentiated axis");

  const auto weights = make_weights(coordinate);

  /*
   * qqはここまでに次元、C連続性、dtype、byteorderを検査済みである。
   * reinterpret_borrowはデータを変換・コピーせず、同じNumPy配列を
   * 要素型Tのpy::array_tとして借り直す。これにより型付きのT*を取得できる。
   */
  auto input = py::reinterpret_borrow<py::array_t<T>>(qq);

  // 入力と同じshape・dtypeの出力を新規確保する。入力データは複製しない。
  py::array_t<T> output({nx, ny, nz});

  const T *input_data = input.data();
  T *output_data = output.mutable_data();

  /*
   * C連続配列ではq[i,j,k]の一次元添字は
   * i*(ny*nz) + j*nz + kとなる。したがって各方向へ1格子進むstrideは
   * x方向がny*nz、y方向がnz、z方向が1である。
   */
  const std::size_t stride_i =
      static_cast<std::size_t>(ny) * static_cast<std::size_t>(nz);
  const std::size_t stride_j = static_cast<std::size_t>(nz);

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
        const std::size_t base =
            static_cast<std::size_t>(i) * stride_i +
            static_cast<std::size_t>(j) * stride_j;

        if constexpr (Axis == 0)
        {
          if (i < margin || i >= nx - margin)
          {
            for (py::ssize_t k = 0; k < nz; ++k)
              output_data[base + static_cast<std::size_t>(k)] = T(0);
            continue;
          }

          const Weights &weight = weights[static_cast<std::size_t>(i)];
          for (py::ssize_t k = 0; k < nz; ++k)
          {
            const std::size_t center = base + static_cast<std::size_t>(k);
            output_data[center] = evaluate(input_data, center, stride_i, weight);
          }
        }
        else if constexpr (Axis == 1)
        {
          if (j < margin || j >= ny - margin)
          {
            for (py::ssize_t k = 0; k < nz; ++k)
              output_data[base + static_cast<std::size_t>(k)] = T(0);
            continue;
          }

          const Weights &weight = weights[static_cast<std::size_t>(j)];
          for (py::ssize_t k = 0; k < nz; ++k)
          {
            const std::size_t center = base + static_cast<std::size_t>(k);
            output_data[center] = evaluate(input_data, center, stride_j, weight);
          }
        }
        else
        {
          output_data[base] = T(0);
          output_data[base + 1] = T(0);
          output_data[base + static_cast<std::size_t>(nz - 2)] = T(0);
          output_data[base + static_cast<std::size_t>(nz - 1)] = T(0);

          for (py::ssize_t k = margin; k < nz - margin; ++k)
          {
            const std::size_t center = base + static_cast<std::size_t>(k);
            output_data[center] = evaluate(
                input_data,
                center,
                1,
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
  auto coordinate_double =
      py::array_t<double, py::array::c_style | py::array::forcecast>::ensure(
          coordinate);
  if (!coordinate_double)
    throw py::type_error("coordinate must be convertible to float64");

  if (is_native_float_dtype<float>(qq.dtype()))
    return derivative_3d<float, Axis>(coordinate_double, qq);
  if (is_native_float_dtype<double>(qq.dtype()))
    return derivative_3d<double, Axis>(coordinate_double, qq);

  throw py::type_error(
      "qq must have native-endian dtype float32 or float64 (" +
      dtype_description(qq.dtype()) + ")");
}
} // namespace derivative_detail
