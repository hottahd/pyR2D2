#pragma once

#include <cstddef>
#include <stdexcept>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <cassert>

namespace py = pybind11;

template <typename T>
struct ViewArray
{
  T *data;
  size_t i_size, j_size, k_size, size;
  size_t ndim;
  inline T &operator()(size_t i)
  {
    return data[i];
  }
  inline const T &operator()(size_t i) const
  {
    return data[i];
  }
  inline T &operator()(size_t i, size_t j)
  {
#ifndef NDEBUG
    assert((ndim >= 2) && "Invalid number of dimensions for 2D access");
#endif
    return data[i * j_size + j];
  }
  inline const T &operator()(size_t i, size_t j) const
  {
#ifndef NDEBUG
    assert((ndim >= 2) && "Invalid number of dimensions for 2D access");
#endif
    return data[i * j_size + j];
  }

  inline T &operator()(size_t i, size_t j, size_t k)
  {
#ifndef NDEBUG
    assert((ndim >= 3) && "Invalid number of dimensions for 3D access");
#endif
    return data[(i * j_size + j) * k_size + k];
  }
  inline const T &operator()(size_t i, size_t j, size_t k) const
  {
#ifndef NDEBUG
    assert((ndim >= 3) && "Invalid number of dimensions for 3D access");
#endif
    return data[(i * j_size + j) * k_size + k];
  }
};

template <typename T>
ViewArray<T> view_array(const py::array_t<T, py::array::c_style | py::array::forcecast> &arr)
{
  auto r = arr.request();

  if (r.ndim < 1 || r.ndim > 3)
    throw std::runtime_error("only 1D, 2D, and 3D arrays are supported");

  ViewArray<T> v;
  v.ndim = static_cast<size_t>(r.ndim);
  v.data = static_cast<T *>(r.ptr);
  v.i_size = static_cast<size_t>(r.shape[0]);

  v.j_size = (r.ndim > 1) ? static_cast<size_t>(r.shape[1]) : 1;
  v.k_size = (r.ndim > 2) ? static_cast<size_t>(r.shape[2]) : 1;
  v.size = static_cast<size_t>(r.size);
  return v;
}