#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <string>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include "numpy_util.hpp"

namespace py = pybind11;

namespace geometry_convert_detail
{
using numpy_detail::CoordinateArray;

/*
 * 現行Fortran版と同じ一様格子を前提とし、格子間隔を検査して返す。
 *
 * spherical2cartesianでは最後から2点目を補間領域上端として使うため、
 * 各球面座標には最低3点を要求する。座標は有限・狭義単調増加でなければ
 * ならない。
 */
inline double uniform_spacing(
    const CoordinateArray &coordinate,
    const char *name)
{
  const py::ssize_t size = coordinate.shape(0);
  if (size < 3)
    throw py::value_error(std::string(name) + " must contain at least 3 points");

  const double *values = coordinate.data();
  const double spacing = values[1] - values[0];
  if (!std::isfinite(values[0]) ||
      !std::isfinite(values[1]) ||
      spacing <= 0.0)
    throw py::value_error(
        std::string(name) + " must be finite and strictly increasing");

  for (py::ssize_t i = 2; i < size; ++i)
  {
    const double current_spacing = values[i] - values[i - 1];
    const double scale = std::max(
        {1.0, std::abs(spacing), std::abs(values[i])});
    // float32座標をfloat64へ変換した場合の丸めも許容する。
    const double tolerance = 1.0e-6 * scale;
    if (!std::isfinite(values[i]) ||
        current_spacing <= 0.0 ||
        std::abs(current_spacing - spacing) > tolerance)
      throw py::value_error(
          std::string(name) +
          " must be finite, strictly increasing, and uniformly spaced");
  }
  return spacing;
}

/*
 * Fortran版と同じく、[-maximum_radius, maximum_radius]を結ぶ一様な
 * 直交座標を逐次加算で生成する。
 */
inline py::array_t<double> make_cartesian_coordinate(
    const py::ssize_t size,
    const double maximum_radius,
    const char *name)
{
  if (size < 2)
    throw py::value_error(std::string(name) + " size must be at least 2");

  py::array_t<double> coordinate(size);
  double *values = coordinate.mutable_data();
  const double spacing =
      (2.0 * maximum_radius) / static_cast<double>(size - 1);
  values[0] = -maximum_radius;
  for (py::ssize_t i = 1; i < size; ++i)
    values[i] = values[i - 1] + spacing;
  return coordinate;
}

/*
 * 一様格子内の値について、下側添字と上側点の線形補間重みを求める。
 *
 * 呼び出し前に値が補間可能領域内であることを検査する。
 */
inline void uniform_bracket(
    const double value,
    const double origin,
    const double spacing,
    std::size_t &lower,
    double &upper_weight)
{
  const double position = (value - origin) / spacing;
  lower = static_cast<std::size_t>(std::floor(position));
  upper_weight = position - static_cast<double>(lower);
}

/*
 * 球面一様格子上の3次元配列を、一様な直交格子へ三線形補間する本体。
 *
 * Tはfloatまたはdoubleで、入力qqsはコピーしない。qqcは入力と同じdtypeの
 * C連続配列として新規確保する。座標、補間重み、8点の積和はdoubleで扱う。
 *
 * 半径と余緯度の有効領域は互換性のため現行Fortran版と同じく、
 * source[0]より大きくsource[-2]より小さい範囲とする。経度は2π周期として
 * 扱い、出力z格子数ではなく球面経度格子数に基づいて折り返す。
 */
template <typename T>
py::tuple spherical_to_cartesian_3d(
    const CoordinateArray &rr,
    const CoordinateArray &th,
    const CoordinateArray &ph,
    const py::array &qqs,
    const py::ssize_t nxc,
    const py::ssize_t nyc,
    const py::ssize_t nzc)
{
  auto input = numpy_detail::borrow_contiguous_3d<T>(qqs, "qqs");

  const py::ssize_t nr = input.shape(0);
  const py::ssize_t nth = input.shape(1);
  const py::ssize_t nph = input.shape(2);
  if (rr.shape(0) != nr || th.shape(0) != nth || ph.shape(0) != nph)
    throw py::value_error(
        "spherical coordinate lengths must match the corresponding qqs axes");

  const double dr = uniform_spacing(rr, "rr");
  const double dth = uniform_spacing(th, "th");
  const double dph = uniform_spacing(ph, "ph");

  const double *radius_values = rr.data();
  const double *theta_values = th.data();
  const double *phi_values = ph.data();
  const double phi_period = phi_values[nph - 1] - phi_values[0];
  const double period_tolerance =
      1.0e-6 * std::max(1.0, 2.0 * std::numbers::pi);
  if (std::abs(phi_period - 2.0 * std::numbers::pi) > period_tolerance)
    throw py::value_error("ph must span one complete 2*pi period");

  const double maximum_radius = radius_values[nr - 1];
  auto xc = make_cartesian_coordinate(nxc, maximum_radius, "xc");
  auto yc = make_cartesian_coordinate(nyc, maximum_radius, "yc");
  auto zc = make_cartesian_coordinate(nzc, maximum_radius, "zc");
  py::array_t<T> output({nxc, nyc, nzc});

  const double *x_values = xc.data();
  const double *y_values = yc.data();
  const double *z_values = zc.data();
  const T *input_data = input.data();
  T *output_data = output.mutable_data();

  const std::size_t input_stride_r =
      numpy_detail::element_stride(input, 0);
  const std::size_t input_stride_th =
      numpy_detail::element_stride(input, 1);
  const std::size_t input_stride_ph =
      numpy_detail::element_stride(input, 2);
  const std::size_t output_stride_x =
      static_cast<std::size_t>(nyc) * static_cast<std::size_t>(nzc);
  const std::size_t output_stride_y = static_cast<std::size_t>(nzc);

  {
    py::gil_scoped_release release;

#pragma omp parallel for collapse(2) schedule(static)
    for (py::ssize_t i = 0; i < nxc; ++i)
    {
      for (py::ssize_t j = 0; j < nyc; ++j)
      {
        const double x = x_values[i];
        const double y = y_values[j];
        const double cylindrical_radius = std::hypot(x, y);
        const std::size_t output_base =
            static_cast<std::size_t>(i) * output_stride_x +
            static_cast<std::size_t>(j) * output_stride_y;

        for (py::ssize_t k = 0; k < nzc; ++k)
        {
          const double z = z_values[k];
          const double radius = std::hypot(cylindrical_radius, z);
          const double theta = std::atan2(cylindrical_radius, z);
          const std::size_t output_index =
              output_base + static_cast<std::size_t>(k);

          if (radius <= radius_values[0] ||
              radius >= radius_values[nr - 2] ||
              theta <= theta_values[0] ||
              theta >= theta_values[nth - 2])
          {
            output_data[output_index] = T(0);
            continue;
          }

          double phi = std::atan2(y, x);
          phi = std::fmod(phi - phi_values[0], phi_period);
          if (phi < 0.0)
            phi += phi_period;
          phi += phi_values[0];

          std::size_t ir0;
          std::size_t ith0;
          std::size_t iph0;
          double wr1;
          double wth1;
          double wph1;
          uniform_bracket(radius, radius_values[0], dr, ir0, wr1);
          uniform_bracket(theta, theta_values[0], dth, ith0, wth1);
          uniform_bracket(phi, phi_values[0], dph, iph0, wph1);

          const std::size_t ir1 = ir0 + 1;
          const std::size_t ith1 = ith0 + 1;
          const std::size_t iph1 = iph0 + 1;
          const double wr0 = 1.0 - wr1;
          const double wth0 = 1.0 - wth1;
          const double wph0 = 1.0 - wph1;

          const auto index = [input_stride_r, input_stride_th, input_stride_ph](
                                 const std::size_t ir,
                                 const std::size_t ith,
                                 const std::size_t iph)
          {
            return ir * input_stride_r +
                   ith * input_stride_th +
                   iph * input_stride_ph;
          };

          const double value =
              static_cast<double>(input_data[index(ir0, ith0, iph0)]) *
                  wr0 * wth0 * wph0 +
              static_cast<double>(input_data[index(ir1, ith0, iph0)]) *
                  wr1 * wth0 * wph0 +
              static_cast<double>(input_data[index(ir0, ith1, iph0)]) *
                  wr0 * wth1 * wph0 +
              static_cast<double>(input_data[index(ir1, ith1, iph0)]) *
                  wr1 * wth1 * wph0 +
              static_cast<double>(input_data[index(ir0, ith0, iph1)]) *
                  wr0 * wth0 * wph1 +
              static_cast<double>(input_data[index(ir1, ith0, iph1)]) *
                  wr1 * wth0 * wph1 +
              static_cast<double>(input_data[index(ir0, ith1, iph1)]) *
                  wr0 * wth1 * wph1 +
              static_cast<double>(input_data[index(ir1, ith1, iph1)]) *
                  wr1 * wth1 * wph1;
          output_data[output_index] = static_cast<T>(value);
        }
      }
    }
  }

  return py::make_tuple(output, xc, yc, zc);
}

/*
 * Python座標を共通形式へ変換し、qqsのdtypeに対応するテンプレート実体を選ぶ。
 */
inline py::tuple spherical_to_cartesian_dispatch(
    const py::array &rr,
    const py::array &th,
    const py::array &ph,
    const py::array &qqs,
    const py::ssize_t nxc,
    const py::ssize_t nyc,
    const py::ssize_t nzc)
{
  const auto rr_double = numpy_detail::ensure_coordinate_1d(rr, "rr");
  const auto th_double = numpy_detail::ensure_coordinate_1d(th, "th");
  const auto ph_double = numpy_detail::ensure_coordinate_1d(ph, "ph");

  if (numpy_detail::is_native_float_dtype<float>(qqs.dtype()))
    return spherical_to_cartesian_3d<float>(
        rr_double, th_double, ph_double, qqs, nxc, nyc, nzc);
  if (numpy_detail::is_native_float_dtype<double>(qqs.dtype()))
    return spherical_to_cartesian_3d<double>(
        rr_double, th_double, ph_double, qqs, nxc, nyc, nzc);

  throw py::type_error(
      "qqs must have native-endian dtype float32 or float64 (" +
      numpy_detail::dtype_description(qqs.dtype()) + ")");
}
} // namespace geometry_convert_detail
