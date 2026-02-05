#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <string>
#include "eos.hpp"
#include "rte.hpp"

namespace py = pybind11;

PYBIND11_MODULE(cpp_util, m)
{
    // clang-format off
  py::class_<EOS>(m, "EOS", R"doc(
  Equation of State (EOS) class for evaluating thermodynamic quantities.

  Parameters
  ----------
  log_ro_np : numpy.ndarray
      1D array of logarithmic density values.
  se_np : numpy.ndarray
      1D array of specific entropy values.
  log_pr_np : numpy.ndarray
      2D array of logarithmic pressure values.
  log_en_np : numpy.ndarray
      2D array of logarithmic internal energy values.
  log_te_np : numpy.ndarray
      2D array of logarithmic temperature values.
  log_op_np : numpy.ndarray
      2D array of logarithmic opacity values.      
    )doc")
      .def(py::init<
           py::array_t<double>,
           py::array_t<double>,
           py::array_t<double>,
           py::array_t<double>,
           py::array_t<double>,
           py::array_t<double>>())
      .def("eval", py::overload_cast<double, double, const std::string &>(&EOS::eval, py::const_), R"doc(
      Evaluate the EOS for given density and specific entropy scalar values.

      Parameters
      ----------
      ro_val : float
          Density value.
      se_val : float
          Specific entropy value.

      var_name : str
          Name of the variable to evaluate. Options are:
          - "pr": Pressure
          - "en": Internal Energy
          - "te": Temperature
          - "op": Opacity
      Returns
      -------
      float
          Evaluated EOS variable.
      )doc")
      .def("eval", py::overload_cast<
        const py::array_t<double, py::array::c_style | py::array::forcecast>&, 
        const py::array_t<double, py::array::c_style | py::array::forcecast>&, 
        const std::string &
        >(&EOS::eval, py::const_), R"doc(
      Evaluate the EOS for given density and specific entropy array values.

      Parameters
      ----------
      ro_val_np : numpy.ndarray
          Array of density values.
      se_val_np : numpy.ndarray
          Array of specific entropy values.
      var_name : str
          Name of the variable to evaluate. Options are:
          - "pr": Pressure
          - "en": Internal Energy
          - "te": Temperature
          - "op": Opacity
      Returns
      -------
      numpy.ndarray
          Array of evaluated EOS variables.
      )doc");
    // clang-format on

    m.def(
        "vertical_upward_rte",
        &vertical_upward_rte,
        R"doc(
    Solve the vertical upward Radiative Transfer Equation (RTE) using the Feautrier method
    for given density, specific entropy, and spatial grids.
    Parameters
    ----------
    ro_np : numpy.ndarray
        3D array of density values.
    se_np : numpy.ndarray
        3D array of specific entropy values.
    x_np : numpy.ndarray
        1D array of spatial grid points.
    eos : EOS
        An instance of the EOS class for thermodynamic evaluations.
    Returns
    -------
    numpy.ndarray
        2D array of radiative transfer results at the top boundary.
    )doc",
        py::arg("ro"),
        py::arg("se"),
        py::arg("x"),
        py::arg("eos"));
}