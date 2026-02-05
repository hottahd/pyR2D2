#pragma once

#include <cmath>
#include "view_array.hpp"

// variable id
static int var_id(const std::string &var_name)
{
  if (var_name == "pr")
    return 0;
  if (var_name == "en")
    return 1;
  if (var_name == "te")
    return 2;
  if (var_name == "op")
    return 3;
  throw std::runtime_error("Unknown variable name: " + var_name);
}

struct EOS
{
  py::array_t<double> log_ro_np;
  py::array_t<double> se_np;

  std::array<py::array_t<double>, 4> log_eos_np;

  size_t i_size;
  size_t j_size;

  double dlog_ro, dse;

  EOS(
      py::array_t<double> log_ro_np_,
      py::array_t<double> se_np_,
      py::array_t<double> log_pr_np_,
      py::array_t<double> log_en_np_,
      py::array_t<double> log_te_np_,
      py::array_t<double> log_op_np_)
      : log_ro_np(log_ro_np_), se_np(se_np_),
        log_eos_np{log_pr_np_, log_en_np_, log_te_np_, log_op_np_}
  {
    auto log_ro = view_array<double>(log_ro_np);
    auto se = view_array<double>(se_np);

    i_size = log_ro.i_size;
    j_size = se.i_size;

    this->dlog_ro = log_ro(1) - log_ro(0);
    this->dse = se(1) - se(0);
  };

  double inline _eval_core(double ro_val, double se_val, const ViewArray<double> &log_ro, const ViewArray<double> &se, const ViewArray<double> &log_eos) const
  {
    double log_ro_val = std::log(ro_val);

    //   find indices
    // size_t i = static_cast<size_t>((log_ro_val - log_ro(0)) / dlog_ro);
    // size_t j = static_cast<size_t>((se_val - se(0)) / dse);

    int ii = (int)std::floor((log_ro_val - log_ro(0)) / dlog_ro);
    int jj = (int)std::floor((se_val - se(0)) / dse);

    if (ii < 0)
      ii = 0;
    if (jj < 0)
      jj = 0;
    if (ii > (int)i_size - 2)
      ii = (int)i_size - 2;
    if (jj > (int)j_size - 2)
      jj = (int)j_size - 2;

    size_t i = (size_t)ii;
    size_t j = (size_t)jj;

    double dlog_ro_val = (log_ro_val - log_ro(i));
    double dse_val = (se_val - se(j));

    // clang-format off
    double qq = std::exp((
      + log_eos(i  , j  ) * (dlog_ro - dlog_ro_val) * (dse - dse_val)
      + log_eos(i+1, j  ) * (          dlog_ro_val) * (dse - dse_val)
      + log_eos(i  , j+1) * (dlog_ro - dlog_ro_val) * (      dse_val)
      + log_eos(i+1, j+1) * (          dlog_ro_val) * (      dse_val)
    )/dlog_ro /dse);
    // clang-format on

    return qq;
  };

  double eval(double ro_val, double se_val, const std::string &var_name) const
  {
    auto log_ro = view_array<double>(log_ro_np);
    auto se = view_array<double>(se_np);
    auto log_eos = view_array<double>(log_eos_np[var_id(var_name)]);

    return _eval_core(ro_val, se_val, log_ro, se, log_eos);
  };

  py::array_t<double> eval(
      const py::array_t<double, py::array::c_style | py::array::forcecast> &ro_val_np,
      const py::array_t<double, py::array::c_style | py::array::forcecast> &se_val_np,
      const std::string &var_name) const
  {
    auto log_ro = view_array<double>(log_ro_np);
    auto se = view_array<double>(se_np);
    auto log_eos = view_array<double>(log_eos_np[var_id(var_name)]);

    auto ro_val = view_array<double>(ro_val_np);
    auto se_val = view_array<double>(se_val_np);

    py::array_t<double> qq_np;
    if (ro_val.ndim == 1)
      qq_np = py::array_t<double>(ro_val.i_size);
    else if (ro_val.ndim == 2)
      qq_np = py::array_t<double>({ro_val.i_size, ro_val.j_size});
    else if (ro_val.ndim == 3)
      qq_np = py::array_t<double>({ro_val.i_size, ro_val.j_size, ro_val.k_size});
    else
      throw std::runtime_error("only 1D, 2D, and 3D arrays are supported");
    auto qq = view_array<double>(qq_np);

    for (size_t idx = 0; idx < ro_val.size; idx++)
    {
      qq(idx) = _eval_core(ro_val(idx), se_val(idx), log_ro, se, log_eos);
    }

    return qq_np;
  };
};
