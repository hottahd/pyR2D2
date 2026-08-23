#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <string>
#include "eos.hpp"
#include "rte.hpp"
#include "yin_yang_convert.hpp"
#include "field_line.hpp"
#include "derivative.hpp"
#include "geometry_convert.hpp"
#include "interpolation.hpp"

namespace py = pybind11;

PYBIND11_MODULE(cpp_util, m)
{
    //
    // Fourth-order derivatives on nonuniform grids
    //
    // clang-format off
    m.def(
        "d_x",
        &derivative_detail::derivative_dispatch<0>,
        py::arg("x"),
        py::arg("qq"),
        R"doc(
        Differentiate a three-dimensional array along the x axis.

        A four-point centered formula is evaluated on a possibly nonuniform
        coordinate grid. The input array is not copied.

        Parameters
        ----------
        x : numpy.ndarray
            One-dimensional x-coordinate array. Its length must equal
            ``qq.shape[0]``. It is converted to C-contiguous float64 when
            necessary.
        qq : numpy.ndarray
            Three-dimensional, C- or F-contiguous, native-endian array with dtype
            float32 or float64 and shape ``(nx, ny, nz)``.

        Returns
        -------
        numpy.ndarray
            C-contiguous derivative array with the same shape and dtype as
            ``qq``.

        Raises
        ------
        TypeError
            If ``qq`` is not native-endian float32 or float64.
        ValueError
            If the array shape, memory layout, coordinate length, or
            coordinate values are invalid.

        Notes
        -----
        The first and last two x planes are set to zero. Coefficients and
        intermediate sums are evaluated in float64. The main loop releases
        the Python GIL and is parallelized with OpenMP.
        )doc");
    m.def(
        "d_y",
        &derivative_detail::derivative_dispatch<1>,
        py::arg("y"),
        py::arg("qq"),
        R"doc(
        Differentiate a three-dimensional array along the y axis.

        A four-point centered formula is evaluated on a possibly nonuniform
        coordinate grid. The input array is not copied.

        Parameters
        ----------
        y : numpy.ndarray
            One-dimensional y-coordinate array. Its length must equal
            ``qq.shape[1]``. It is converted to C-contiguous float64 when
            necessary.
        qq : numpy.ndarray
            Three-dimensional, C- or F-contiguous, native-endian array with dtype
            float32 or float64 and shape ``(nx, ny, nz)``.

        Returns
        -------
        numpy.ndarray
            C-contiguous derivative array with the same shape and dtype as
            ``qq``.

        Raises
        ------
        TypeError
            If ``qq`` is not native-endian float32 or float64.
        ValueError
            If the array shape, memory layout, coordinate length, or
            coordinate values are invalid.

        Notes
        -----
        The first and last two y planes are set to zero. Coefficients and
        intermediate sums are evaluated in float64. The main loop releases
        the Python GIL and is parallelized with OpenMP.
        )doc");
    m.def(
        "d_z",
        &derivative_detail::derivative_dispatch<2>,
        py::arg("z"),
        py::arg("qq"),
        R"doc(
        Differentiate a three-dimensional array along the z axis.

        A four-point centered formula is evaluated on a possibly nonuniform
        coordinate grid. The input array is not copied.

        Parameters
        ----------
        z : numpy.ndarray
            One-dimensional z-coordinate array. Its length must equal
            ``qq.shape[2]``. It is converted to C-contiguous float64 when
            necessary.
        qq : numpy.ndarray
            Three-dimensional, C- or F-contiguous, native-endian array with dtype
            float32 or float64 and shape ``(nx, ny, nz)``.

        Returns
        -------
        numpy.ndarray
            C-contiguous derivative array with the same shape and dtype as
            ``qq``.

        Raises
        ------
        TypeError
            If ``qq`` is not native-endian float32 or float64.
        ValueError
            If the array shape, memory layout, coordinate length, or
            coordinate values are invalid.

        Notes
        -----
        The first and last two z planes are set to zero. Coefficients and
        intermediate sums are evaluated in float64. The main loop releases
        the Python GIL and is parallelized with OpenMP.
        )doc");
    m.def(
        "interp",
        &interpolation_detail::interpolation_dispatch,
        py::arg("x"),
        py::arg("y"),
        py::arg("z"),
        py::arg("xu"),
        py::arg("yu"),
        py::arg("zu"),
        py::arg("qq"),
        R"doc(
        Trilinearly interpolate a three-dimensional array.

        The source grid may be nonuniform but each source coordinate must be
        finite and strictly increasing. The input array is not copied.

        Parameters
        ----------
        x, y, z : numpy.ndarray
            One-dimensional source coordinates. Their lengths must match the
            corresponding axes of ``qq``.
        xu, yu, zu : numpy.ndarray
            One-dimensional target coordinates.
        qq : numpy.ndarray
            Three-dimensional, C- or F-contiguous, native-endian array with dtype
            float32 or float64 and shape ``(len(x), len(y), len(z))``.

        Returns
        -------
        numpy.ndarray
            C-contiguous array with shape ``(len(xu), len(yu), len(zu))`` and
            the same dtype as ``qq``.

        Raises
        ------
        TypeError
            If ``qq`` is not native-endian float32 or float64.
        ValueError
            If an input shape, memory layout, or coordinate is invalid.

        Notes
        -----
        Target points on the source-domain boundary or outside the source
        domain are set to zero, matching the existing Fortran implementation.
        Weights and intermediate sums are evaluated in float64. The main loop
        releases the Python GIL and is parallelized with OpenMP.
        )doc");
    m.def(
        "spherical2cartesian",
        &geometry_convert_detail::spherical_to_cartesian_dispatch,
        py::arg("rr"),
        py::arg("th"),
        py::arg("ph"),
        py::arg("qqs"),
        py::arg("ixc"),
        py::arg("jxc"),
        py::arg("kxc"),
        R"doc(
        Interpolate a scalar field from a spherical grid to a Cartesian grid.

        Parameters
        ----------
        rr, th, ph : numpy.ndarray
            One-dimensional, uniformly spaced spherical coordinates. ``ph``
            must span one complete 2*pi period.
        qqs : numpy.ndarray
            Three-dimensional, C- or F-contiguous, native-endian float32 or float64
            array with shape ``(len(rr), len(th), len(ph))``.
        ixc, jxc, kxc : int
            Numbers of points in the Cartesian x, y, and z coordinates. Each
            must be at least two.

        Returns
        -------
        qqc : numpy.ndarray
            C-contiguous Cartesian scalar field with the same dtype as
            ``qqs`` and shape ``(ixc, jxc, kxc)``.
        xc, yc, zc : numpy.ndarray
            Float64 Cartesian coordinates spanning ``[-max(rr), max(rr)]``.

        Notes
        -----
        The radial and colatitude interpolation domain matches the existing
        Fortran implementation: values must lie strictly between the first
        and second-to-last source coordinates. Other output points are zero.
        Longitude is treated as periodic. The input field is not copied,
        calculations release the Python GIL, and the output loop is
        parallelized with OpenMP.
        )doc");
    // clang-format on


    //
    // EOS
    //
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
           py::array_t<float>,
           py::array_t<float>,
           py::array_t<float>,
           py::array_t<float>,
           py::array_t<float>,
           py::array_t<float>>())
      .def("eval", py::overload_cast<float, float, const std::string &>(&EOS::eval, py::const_), R"doc(
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
        const py::array_t<float, py::array::c_style | py::array::forcecast>&, 
        const py::array_t<float, py::array::c_style | py::array::forcecast>&, 
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

      Warning
      -------
      This reads the table unconditionally, which is not what R2D2 does. The
      code switches to the linearized EOS whenever the relative amplitude
      ct = max(|ro1|/ro0, |se1|/se0) stays below ct0 = 3e-3, and in a
      deep-convection run that is everywhere. Use pyR2D2.util.eos_switch to
      reproduce what the code actually did.
      )doc");
    // clang-format on

    //
    // Radiative transfer
    //

    m.def(
        "eval_tau",
        &eval_tau,
        R"doc(
    Evaluate the optical depth for given density, specific entropy, and spatial grids.
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
    i_bot : int, optional
        The bottom index for the evaluation range (default is 0).
    Returns
    -------
    numpy.ndarray
        2D array of radiative transfer results at the top boundary.
    )doc",
        py::arg("ro"),
        py::arg("se"),
        py::arg("x"),
        py::arg("eos"),
        py::arg("i_bot") = 0);

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

    //
    // Yin-Yang grid conversion
    //

    py::class_<YinYang>(m, "YinYang", R"doc(
    Class for converting data between Yin-Yang grid and regular spherical grid.
    Parameters
    ----------
    th_yy_np : numpy.ndarray
        1D array of Yin-Yang grid theta values.
    ph_yy_np : numpy.ndarray
        1D array of Yin-Yang grid phi values.
    th_np : numpy.ndarray
        1D array of regular grid theta values.
    ph_np : numpy.ndarray
        1D array of regular grid phi values.
     )doc")
        .def(py::init<
             py::array_t<double, py::array::c_style | py::array::forcecast>,
             py::array_t<double, py::array::c_style | py::array::forcecast>,
             py::array_t<double, py::array::c_style | py::array::forcecast>,
             py::array_t<double, py::array::c_style | py::array::forcecast>>())
        .def(
            "convert_scalar",
            &YinYang::convert_scalar<float>,
            R"doc(
    Convert scalar values from Yin-Yang grid to regular grid.
    
    Parameters
    ----------
    qq_yin : numpy.ndarray (dtype=float32)
        2D or 3D array of Yin grid values. If 3D, the last index indicates the variable index.
    qq_yan : numpy.ndarray (dtype=float32)
        2D or 3D array of Yang grid values. If 3D, the last index indicates the variable index.

    Returns
    -------
    None
    )doc",
            py::arg("qq_yin"),
            py::arg("qq_yan"));

    //
    // Field line tracing
    //

    m.def(
        "trace_field_line",
        &trace_field_line,
        R"doc(
    Trace a field line through a 3D magnetic field.
    Parameters
    ----------
    x, y, z : numpy.ndarray
        1D coordinate arrays. Must be uniform (checked, with an error rather
        than silent corruption). Read as float64; positions are integrated in
        double as offsets from x[0], y[0], z[0], so a large absolute offset no
        longer eats into the precision.
    bx, by, bz : numpy.ndarray
        3D magnetic field components.
    fields : dict
        Extra 3D arrays to interpolate along the line.
    x0, y0, z0 : double
        Starting point, in the same absolute coordinates as x, y, z.
    ds : double
        Arc-length step.
    n_steps : int
        Maximum number of points.
    sign : double, optional
        +1 traces along B, -1 against it. Default +1.
    method : str, optional
        'implicit_midpoint' (default) is time symmetric, so tracing back with
        -ds retraces the same points to tolerance. 'explicit_midpoint' is the
        pre-2026-08-23 behaviour, kept for comparison.
    max_iter : int, optional
        Fixed-point iterations allowed per step for the implicit scheme.
    tol : double, optional
        Convergence tolerance, relative to ds.

    Returns
    -------
    dict
        'x', 'y', 'z' : float64 positions, truncated at n_valid
        each key of `fields` : float32 values along the line
        'n_valid' : number of points produced
        'status'  : 0 ran out of steps, 1 left the domain, 2 |B| vanished,
                    3 iteration did not converge
    )doc",
        py::arg("x"),
        py::arg("y"),
        py::arg("z"),
        py::arg("bx"),
        py::arg("by"),
        py::arg("bz"),
        py::arg("fields"),
        py::arg("x0"),
        py::arg("y0"),
        py::arg("z0"),
        py::arg("ds"),
        py::arg("n_steps"),
        py::arg("sign") = 1.0,
        py::arg("method") = "implicit_midpoint",
        py::arg("max_iter") = 20,
        py::arg("tol") = 1e-10);

    m.def(
        "trace_field_lines",
        &trace_field_lines,
        R"doc(
    Trace many field lines at once (OpenMP parallel).

    Same physics as :func:`trace_field_line`, but the seeds are arrays and the
    output is rectangular, which avoids one Python call per line.

    Parameters
    ----------
    x, y, z : numpy.ndarray
        1D coordinate arrays. Must be uniform (checked).
    bx, by, bz : numpy.ndarray
        3D magnetic field components.
    fields : dict
        Extra 3D arrays to interpolate along each line.
    x0, y0, z0 : numpy.ndarray
        1D arrays of seed coordinates, one entry per line.
    ds : double
        Arc-length step.
    n_steps : int
        Maximum number of points per line.
    sign : double, optional
        +1 traces along B, -1 against it. Default +1.
    method : str, optional
        'implicit_midpoint' (default, time symmetric) or 'explicit_midpoint'
        (the pre-2026-08-23 behaviour, kept for comparison).
    max_iter : int, optional
        Fixed-point iterations allowed per step for the implicit scheme.
    tol : double, optional
        Convergence tolerance, relative to ds.

    Returns
    -------
    dict
        'x', 'y', 'z' : (n_lines, n_steps) float64, NaN past the end of a line
        each key of `fields` : (n_lines, n_steps) float32, NaN past the end
        'n_valid' : (n_lines,) int64, number of populated points
        'status'  : (n_lines,) int32, 0 ran out of steps, 1 left the domain,
                    2 |B| vanished, 3 iteration did not converge
    )doc",
        py::arg("x"),
        py::arg("y"),
        py::arg("z"),
        py::arg("bx"),
        py::arg("by"),
        py::arg("bz"),
        py::arg("fields"),
        py::arg("x0"),
        py::arg("y0"),
        py::arg("z0"),
        py::arg("ds"),
        py::arg("n_steps"),
        py::arg("sign") = 1.0,
        py::arg("method") = "implicit_midpoint",
        py::arg("max_iter") = 20,
        py::arg("tol") = 1e-10);
}
