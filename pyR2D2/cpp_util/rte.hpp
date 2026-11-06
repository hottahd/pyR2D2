#pragma once
#include <vector>
#include <iostream>
#include <cmath>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "const.hpp"
#include "view_array.hpp"
#include "eos.hpp"

namespace rte_detail
{
  // Below these thresholds the general formulas lose all their digits (0/0 at worst), and the series is exact to
  // rounding: the error of the first terms is O(d^2) ~ 1e-12 relative. Cf. rte_cell_crossing.hpp of R2D2plus.
  constexpr double small_difference = 1.0e-6;

  // Optical thickness of a cell where log(alpha) is linear in x
  inline double cell_dtu(double al_end, double al_stt, double log_al_end, double log_al_stt, double dx)
  {
    const double dlog = log_al_end - log_al_stt;
    const double mean = std::abs(dlog) < small_difference ? 0.5 * (al_end + al_stt) : (al_end - al_stt) / dlog;
    return mean * dx;
  }

  // Intensity added in a cell where log(S) is linear in tau
  inline double cell_din(double so_end, double so_stt, double log_so_end, double log_so_stt, double dtu)
  {
    const double attenuation = std::exp(-dtu);
    const double denominator = log_so_end - log_so_stt + dtu;
    if (std::abs(denominator) < small_difference)
    {
      return dtu * so_stt * attenuation * (1.0 + 0.5 * denominator);
    }
    return dtu * (so_end - so_stt * attenuation) / denominator;
  }
} // namespace rte_detail

py::array_t<float> eval_tau(
    const py::array_t<float, py::array::c_style> &ro_np,
    const py::array_t<float, py::array::c_style> &se_np,
    const py::array_t<double, py::array::c_style | py::array::forcecast> &x_np,
    EOS &eos,
    size_t i_bot = 0)
{

  auto ro = view_array<float>(ro_np);

  if (i_bot >= ro.i_size)
  {
    throw std::runtime_error("i_bot must be smaller than ro.i_size");
  }

  auto x = view_array<double>(x_np);

  if (x.i_size != ro.i_size)
  {
    throw std::runtime_error("x_np length must match ro_np's first axis");
  }

  py::array_t<float> tu_np({ro.i_size, ro.j_size, ro.k_size});
  auto tu = view_array<float>(tu_np);

  py::array_t<float> op_np = eos.eval(ro_np, se_np, "op");
  auto op = view_array<float>(op_np);

  {
    // No Python API below: only raw pointers
    py::gil_scoped_release release;

#pragma omp parallel
    {
      std::vector<double> log_al(ro.i_size), al(ro.i_size);

#pragma omp for collapse(2) schedule(static)
      for (size_t j = 0; j < ro.j_size; j++)
      {
        for (size_t k = 0; k < ro.k_size; k++)
        {

          for (size_t i = i_bot; i < ro.i_size; i++)
          {
            al[i] = static_cast<double>(ro(i, j, k)) * static_cast<double>(op(i, j, k));
            log_al[i] = std::log(al[i]);
          }

          double tu_sum = 0.0;
          tu(ro.i_size - 1, j, k) = 0.0f;
          for (size_t i = ro.i_size - 1; i > i_bot; --i)
          {
            tu_sum += rte_detail::cell_dtu(al[i - 1], al[i], log_al[i - 1], log_al[i], std::abs(x(i - 1) - x(i)));
            tu(i - 1, j, k) = static_cast<float>(tu_sum);
          }

          for (size_t i = 0; i < i_bot; i++)
          {
            tu(i, j, k) = tu(i_bot, j, k);
          }
        }
      }
    }
  }

  return tu_np;
}

py::array_t<float> vertical_upward_rte(
    const py::array_t<float, py::array::c_style> &ro_np,
    const py::array_t<float, py::array::c_style> &se_np,
    const py::array_t<double, py::array::c_style | py::array::forcecast> &x_np,
    EOS &eos)
{

  auto ro = view_array<float>(ro_np);
  auto se = view_array<float>(se_np);
  auto x = view_array<double>(x_np);

  if (x.i_size != ro.i_size)
  {
    throw std::runtime_error("x_np length must match ro_np's first axis");
  }
  if (ro.i_size < 2)
  {
    throw std::runtime_error("ro_np's first axis must have at least 2 cells");
  }

  py::array_t<float> rt_np({ro.j_size, ro.k_size});
  auto rt = view_array<float>(rt_np);

  py::array_t<float> op_np = eos.eval(ro_np, se_np, "op");
  py::array_t<float> te_np = eos.eval(ro_np, se_np, "te");

  auto op = view_array<float>(op_np);
  auto te = view_array<float>(te_np);

  std::vector<double> x_mid(ro.i_size);
  for (size_t i = 1; i < ro.i_size; i++)
  {
    x_mid[i] = 0.5 * (x(i - 1) + x(i));
  };

  {
    // No Python API below: only raw pointers
    py::gil_scoped_release release;

#pragma omp parallel
    {
      std::vector<double> log_al_ctr(ro.i_size), log_so_ctr(ro.i_size);
      std::vector<double> log_al_mid(ro.i_size), log_so_mid(ro.i_size);
      std::vector<double> al_mid(ro.i_size), so_mid(ro.i_size);

#pragma omp for collapse(2) schedule(static)
      for (size_t j = 0; j < ro.j_size; j++)
      {
        for (size_t k = 0; k < ro.k_size; k++)
        {
          for (size_t i = 0; i < ro.i_size; i++)
          {
            log_al_ctr[i] = std::log(static_cast<double>(ro(i, j, k)) * static_cast<double>(op(i, j, k)));
            log_so_ctr[i] = std::log(cst::sb * std::pow(static_cast<double>(te(i, j, k)), 4) * cst::pii);
          }

          // value at i - 1/2
          for (size_t i = 1; i < ro.i_size; i++)
          {
            log_al_mid[i] = 0.5 * (log_al_ctr[i - 1] + log_al_ctr[i]);
            log_so_mid[i] = 0.5 * (log_so_ctr[i - 1] + log_so_ctr[i]);
            al_mid[i] = std::exp(log_al_mid[i]);
            so_mid[i] = std::exp(log_so_mid[i]);
          }

          // From the lowest cell face (I = S) to the highest one, i_size - 1 - 1/2
          double intensity = so_mid[1];
          for (size_t i = 1; i < ro.i_size - 1; i++)
          {
            const double dtu = rte_detail::cell_dtu(
                al_mid[i + 1], al_mid[i], log_al_mid[i + 1], log_al_mid[i], x_mid[i + 1] - x_mid[i]);
            intensity = intensity * std::exp(-dtu) + rte_detail::cell_din(so_mid[i + 1], so_mid[i], log_so_mid[i + 1], log_so_mid[i], dtu);
          }
          rt(j, k) = static_cast<float>(intensity);
        }
      }
    }
  }

  return rt_np;
}
