#pragma once
#include <vector>
#include <iostream>
#include <cmath>
#include <omp.h>
#include "const.hpp"
#include "view_array.hpp"
#include "eos.hpp"

py::array_t<float> eval_tau(
    const py::array_t<float, py::array::c_style> &ro_np,
    const py::array_t<float, py::array::c_style> &se_np,
    const py::array_t<float, py::array::c_style> &x_np,
    EOS &eos,
    size_t i_bot = 0)
{

  auto ro = view_array<float>(ro_np);

  if (i_bot >= ro.i_size)
  {
    throw std::runtime_error("i_bot must be smaller than ro.i_size");
  }

  auto x = view_array<float>(x_np);

  py::array_t<float> tu_np({ro.i_size, ro.j_size, ro.k_size});
  auto tu = view_array<float>(tu_np);

  py::array_t<float> op_np = eos.eval(ro_np, se_np, "op");
  auto op = view_array<float>(op_np);

  // py::gil_scoped_release release;

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
          al[i] = ro(i, j, k) * op(i, j, k);
          log_al[i] = std::log(al[i]);
        }

        tu(ro.i_size - 1, j, k) = 0.0;
        for (size_t i = ro.i_size - 1; i > i_bot; --i)
        {
          float log_al_stt = log_al[i];
          float al_stt = al[i];
          float x_stt = x(i);

          float log_al_end = log_al[i - 1];
          float al_end = al[i - 1];
          float x_end = x(i - 1);

          float dtu = (al_end - al_stt) * std::abs(x_end - x_stt) / (log_al_end - log_al_stt);
          tu(i - 1, j, k) = tu(i, j, k) + dtu;
        }

        for (size_t i = 0; i < i_bot; i++)
        {
          tu(i, j, k) = tu(i_bot, j, k);
        }
      }
    }
  }

  return tu_np;
}

py::array_t<float> vertical_upward_rte(
    const py::array_t<float, py::array::c_style> &ro_np,
    const py::array_t<float, py::array::c_style> &se_np,
    const py::array_t<float, py::array::c_style> &x_np,
    EOS &eos)
{

  auto ro = view_array<float>(ro_np);
  auto se = view_array<float>(se_np);
  auto x = view_array<float>(x_np);

  py::array_t<float> rt_np({ro.i_size, ro.j_size, ro.k_size});
  auto rt = view_array<float>(rt_np);

  py::array_t<float> op_np = eos.eval(ro_np, se_np, "op");
  py::array_t<float> te_np = eos.eval(ro_np, se_np, "te");

  auto op = view_array<float>(op_np);
  auto te = view_array<float>(te_np);

  std::vector<float> x_mid(ro.i_size);
  for (size_t i = 1; i < ro.i_size; i++)
  {
    x_mid[i] = 0.5 * (x(i - 1) + x(i));
  };

  // py::gil_scoped_release release;

#pragma omp parallel
  {
    std::vector<float> log_al_ctr(ro.i_size), log_so_ctr(ro.i_size);
    std::vector<float> log_al_mid(ro.i_size), log_so_mid(ro.i_size);
    std::vector<float> al_mid(ro.i_size), so_mid(ro.i_size);
    std::vector<float> rt1d(ro.i_size);

#pragma omp for collapse(2) schedule(static)
    for (size_t j = 0; j < ro.j_size; j++)
    {
      for (size_t k = 0; k < ro.k_size; k++)
      {
        for (size_t i = 0; i < ro.i_size; i++)
        {
          log_al_ctr[i] = std::log(ro(i, j, k) * op(i, j, k));
          log_so_ctr[i] = std::log(cst::sb * std::pow(te(i, j, k), 4) * cst::pii);
        }

        // value at i - 1/2
        for (size_t i = 1; i < ro.i_size; i++)
        {
          log_al_mid[i] = 0.5 * (log_al_ctr[i - 1] + log_al_ctr[i]);
          log_so_mid[i] = 0.5 * (log_so_ctr[i - 1] + log_so_ctr[i]);
          al_mid[i] = std::exp(log_al_mid[i]);
          so_mid[i] = std::exp(log_so_mid[i]);
        }

        rt1d[1] = so_mid[1];
        for (size_t i = 1; i < ro.i_size - 1; i++)
        {
          float log_so_end = log_so_mid[i + 1];
          float log_al_end = log_al_mid[i + 1];
          float so_end = so_mid[i + 1];
          float al_end = al_mid[i + 1];
          float x_end = x_mid[i + 1];

          float log_so_stt = log_so_mid[i];
          float log_al_stt = log_al_mid[i];
          float so_stt = so_mid[i];
          float al_stt = al_mid[i];
          float x_stt = x_mid[i];

          float dtu = (al_end - al_stt) * (x_end - x_stt) / (log_al_end - log_al_stt);
          // tu(i, j, k) = dtu;

          rt1d[i + 1] = rt1d[i] * std::exp(-dtu) + dtu * (so_end - so_stt * std::exp(-dtu)) / (log_so_end - log_so_stt + dtu);
        }
        rt(j, k) = rt1d[ro.i_size - 2];
      }
    }
  }

  return rt_np;
}