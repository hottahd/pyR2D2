#pragma once

#include <bit>
#include <cstddef>
#include <string>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace numpy_detail
{
/*
 * 1次元座標に共通して使用するNumPy配列型。
 *
 * 巨大な3次元物理量とは異なり、座標は十分小さいため、必要に応じた
 * float64変換とC連続化を許可する。
 */
using CoordinateArray =
    py::array_t<double, py::array::c_style | py::array::forcecast>;

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
 * Pythonから渡された座標を、C連続な1次元float64配列として用意する。
 *
 * nameはエラーメッセージに使用し、x、xuなど、どの座標で失敗したかを
 * 呼び出し側が判別できるようにする。
 */
inline CoordinateArray ensure_coordinate_1d(
    const py::array &coordinate,
    const char *name)
{
  auto converted = CoordinateArray::ensure(coordinate);
  if (!converted)
    throw py::type_error(
        std::string(name) + " must be convertible to a float64 NumPy array");
  if (converted.ndim() != 1)
    throw py::value_error(std::string(name) + " must be a 1D array");
  return converted;
}

/*
 * 巨大な3次元NumPy配列を、コピーせずpy::array_t<T>として借りる。
 *
 * reinterpret_borrowを行う前に、次元、CまたはF連続性、dtype、byteorderを検査する。
 * 条件を満たさない配列を暗黙変換しないため、入力と同規模の一時配列は
 * ここでは作られない。
 */
template <typename T>
inline py::array_t<T> borrow_contiguous_3d(
    const py::array &array,
    const char *name)
{
  if (array.ndim() != 3)
    throw py::value_error(std::string(name) + " must be a 3D array");
  const bool c_contiguous = array.flags() & py::array::c_style;
  const bool f_contiguous = array.flags() & py::array::f_style;
  if (!c_contiguous && !f_contiguous)
    throw py::value_error(std::string(name) + " must be C- or F-contiguous");
  if (!is_native_float_dtype<T>(array.dtype()))
    throw py::type_error(
        std::string(name) +
        " must have native-endian dtype matching the selected kernel (" +
        dtype_description(array.dtype()) + ")");

  return py::reinterpret_borrow<py::array_t<T>>(array);
}
/*
 * NumPyのバイト単位strideを、Tの要素数単位へ変換する。
 *
 * C連続・F連続のどちらでも、この値を使えば同じ添字式で入力を読める。
 */
template <typename T>
inline std::size_t element_stride(
    const py::array_t<T> &array,
    const py::ssize_t axis)
{
  return static_cast<std::size_t>(array.strides(axis)) / sizeof(T);
}
} // namespace numpy_detail
