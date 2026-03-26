#pragma once

#include <cstddef>
#include <stdexcept>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <cassert>

namespace py = pybind11;

// Lower-dimensional arrays can be accessed with higher-dimensional operators
// if extra indices are zero:
// 1D: a(i), a(i,0), a(i,0,0)
// 2D: a(i,j), a(i,j,0)
// 3D: a(i,j,k)
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
    if (ndim == 1)
      assert(j == 0 && "For 1D access with (i,j), j must be 0");
#endif
    return data[i * j_size + j];
  }
  inline const T &operator()(size_t i, size_t j) const
  {
#ifndef NDEBUG
    if (ndim == 1)
      assert(j == 0 && "For 1D access with (i,j), j must be 0");
#endif
    return data[i * j_size + j];
  }

  inline T &operator()(size_t i, size_t j, size_t k)
  {
#ifndef NDEBUG
    if (ndim == 1)
    {
      assert(j == 0 && "For 1D access with (i,j,k), j must be 0");
      assert(k == 0 && "For 1D access with (i,j,k), k must be 0");
    }
    if (ndim == 2)
      assert(k == 0 && "For 2D access with (i,j,k), k must be 0");
#endif
    return data[(i * j_size + j) * k_size + k];
  }
  inline const T &operator()(size_t i, size_t j, size_t k) const
  {
#ifndef NDEBUG
    if (ndim == 1)
    {
      assert(j == 0 && "For 1D access with (i,j,k), j must be 0");
      assert(k == 0 && "For 1D access with (i,j,k), k must be 0");
    }
    if (ndim == 2)
      assert(k == 0 && "For 2D access with (i,j,k), k must be 0");
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