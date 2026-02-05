#include <vector>
#include <iostream>
#include <cmath>
#include "view_array.hpp"
#include "eos.hpp"

double pi = 3.14159265358979323846;
double pii = 1.0 / pi;

py::array_t<double> vertical_upward_rte(
    const py::array_t<double, py::array::c_style | py::array::forcecast> &ro_np,
    const py::array_t<double, py::array::c_style | py::array::forcecast> &se_np,
    const py::array_t<double, py::array::c_style | py::array::forcecast> &x_np,
    EOS &eos)
{

  double sb = 5.67e-5; // Stefan-Boltzmann constant in cgs units
  auto ro = view_array<double>(ro_np);
  auto se = view_array<double>(se_np);
  auto x = view_array<double>(x_np);

  py::array_t<double> rt_np({ro.j_size, ro.k_size});
  auto rt = view_array<double>(rt_np);
  py::array_t<double> op_np = eos.eval(ro_np, se_np, "op");
  py::array_t<double> te_np = eos.eval(ro_np, se_np, "te");

  auto op = view_array<double>(op_np);
  auto te = view_array<double>(te_np);

  std::vector<double> x_mid(ro.i_size);
  for (size_t i = 1; i < ro.i_size; i++)
  {
    x_mid[i] = 0.5 * (x(i - 1) + x(i));
  };

  std::vector<double> log_al_ctr(ro.i_size), log_so_ctr(ro.i_size);
  std::vector<double> log_al_mid(ro.i_size), log_so_mid(ro.i_size);
  std::vector<double> al_mid(ro.i_size), so_mid(ro.i_size);
  std::vector<double> rt1d(ro.i_size);

  for (size_t j = 0; j < ro.j_size; j++)
  {
    for (size_t k = 0; k < ro.k_size; k++)
    {
      for (size_t i = 0; i < ro.i_size; i++)
      {
        log_al_ctr[i] = std::log(ro(i, j, k) * op(i, j, k));
        log_so_ctr[i] = std::log(sb * std::pow(te(i, j, k), 4) * pii);
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
        double log_so_end = log_so_mid[i + 1];
        double log_al_end = log_al_mid[i + 1];
        double so_end = so_mid[i + 1];
        double al_end = al_mid[i + 1];
        double x_end = x_mid[i + 1];

        double log_so_stt = log_so_mid[i];
        double log_al_stt = log_al_mid[i];
        double so_stt = so_mid[i];
        double al_stt = al_mid[i];
        double x_stt = x_mid[i];

        double dtu = (al_end - al_stt) * (x_end - x_stt) / (log_al_end - log_al_stt);

        rt1d[i + 1] = rt1d[i] * std::exp(-dtu) + dtu * (so_end - so_stt * std::exp(-dtu)) / (log_so_end - log_so_stt + dtu);
      }
      rt(j, k) = rt1d[ro.i_size - 2];
    }
  }

  return rt_np;
}