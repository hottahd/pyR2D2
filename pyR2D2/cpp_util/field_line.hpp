#pragma once
#include <cmath>
#include <cstdint>
#include <string>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "view_array.hpp"

// Field line tracing.
//
// Public entry points:
//   trace_field_line   -- one line, returns 1D arrays truncated at n_valid
//   trace_field_lines  -- many lines at once, OpenMP parallel, 2D arrays
//
// Design notes (2026-08-23 rewrite)
// --------------------------------
// * Positions are integrated in **double**, as offsets from the grid origin.
//   The previous version integrated in float32 on the caller's absolute
//   coordinates.  With y ~ 2.4e9 cm the float32 spacing is 256 cm, i.e. 0.26 %
//   of a 1 km step, which random-walks to ~0.4 km over 20000 steps -- larger
//   than the measured round-trip error of the scheme itself.  Worse, a caller
//   passing an absolute x ~ 6.9e10 cm got a spacing of 8192 cm (8 % of a step,
//   ~12 km over 20000 steps) with no warning.  Subtracting the origin inside
//   the function removes that trap entirely.
//
// * The default integrator is the **implicit midpoint rule**
//       p_{n+1} = p_n + ds * u((p_n + p_{n+1})/2),      u = B/|B|
//   solved by fixed-point iteration.  It is *time symmetric*: tracing back
//   with -ds from p_{n+1} returns to p_n to iteration tolerance.  The previous
//   explicit midpoint rule is not, so a down-then-up round trip accumulated a
//   scheme error that was indistinguishable from genuine chaotic divergence of
//   the field.  With a reversible scheme the round-trip error becomes a clean
//   diagnostic of the field itself.  The old scheme is still reachable with
//   method="explicit_midpoint" for comparison.
//
// * |B| = 0 is now guarded (the old code divided unconditionally), the grid is
//   checked for uniformity (it was assumed silently), and every trace reports
//   why it stopped.
//
// Interpolation is still trilinear, i.e. C0.  That caps the effective order of
// any integrator at ~2, so moving to RK4 alone would buy little; tricubic is
// the next step if the numbers call for it.

struct FieldData
{
    std::string name;

    py::array_t<float, py::array::c_style | py::array::forcecast> qq_np;
    py::array_t<float> out_np;

    ViewArray<float> qq;
    ViewArray<float> out;
};

// Why a trace stopped.
enum TraceStatus : int32_t
{
    TRACE_OK = 0,             // ran out of n_steps
    TRACE_OUT_OF_DOMAIN = 1,  // the line (or its midpoint) left the grid
    TRACE_ZERO_FIELD = 2,     // |B| vanished, direction undefined
    TRACE_NO_CONVERGE = 3,    // fixed-point iteration hit max_iter
};

// A uniform 1D grid, held in double.
struct Grid1D
{
    double origin;
    double delta;
    double span;  // (n - 2) * delta: last position interpolate3d can read from
    size_t n;
};

inline Grid1D make_grid(const ViewArray<double> &a, const char *name)
{
    if (a.i_size < 2)
    {
        throw std::runtime_error(std::string("trace_field_line: ") + name +
                                 " grid needs at least 2 points");
    }
    const double d0 = double(a(1)) - double(a(0));
    if (!(std::abs(d0) > 0.0))
    {
        throw std::runtime_error(std::string("trace_field_line: ") + name +
                                 " grid spacing is zero");
    }
    // The interpolation assumes a uniform grid.  R2D2 allows a non-uniform x
    // grid, so check rather than break silently.
    const double tol = 1e-3 * std::abs(d0);
    for (size_t i = 1; i + 1 < a.i_size; i++)
    {
        if (std::abs((double(a(i + 1)) - double(a(i))) - d0) > tol)
        {
            throw std::runtime_error(
                std::string("trace_field_line: ") + name +
                " grid is not uniform; this routine assumes a uniform grid. "
                "(If the array is float32 with a large absolute offset -- e.g. an "
                "R2D2 physical x of ~7e10 cm -- its spacing cannot be represented "
                "exactly; pass float64, or coordinates relative to the domain.)");
        }
    }
    Grid1D g;
    g.origin = double(a(0));
    g.delta = d0;
    g.n = a.i_size;
    g.span = double(a.i_size - 2) * d0;
    return g;
}

// Trilinear interpolation at a position given as an offset from the grid
// origin.  Callers must have checked that the offset is inside [0, span].
inline double interpolate3d_offset(const ViewArray<float> &qq,
                                   double ox, double oy, double oz,
                                   const Grid1D &gx, const Grid1D &gy, const Grid1D &gz)
{
    double fx = ox / gx.delta;
    double fy = oy / gy.delta;
    double fz = oz / gz.delta;

    long i = (long)std::floor(fx);
    long j = (long)std::floor(fy);
    long k = (long)std::floor(fz);

    // Guard against the floor landing one past the end at the very edge.
    if (i < 0) i = 0;
    if (j < 0) j = 0;
    if (k < 0) k = 0;
    if ((size_t)i > gx.n - 2) i = (long)gx.n - 2;
    if ((size_t)j > gy.n - 2) j = (long)gy.n - 2;
    if ((size_t)k > gz.n - 2) k = (long)gz.n - 2;

    const double dx1 = fx - double(i), dx2 = 1.0 - dx1;
    const double dy1 = fy - double(j), dy2 = 1.0 - dy1;
    const double dz1 = fz - double(k), dz2 = 1.0 - dz1;

    // clang-format off
    return
    double(qq(i    , j    , k    ))*dx2*dy2*dz2 +
    double(qq(i + 1, j    , k    ))*dx1*dy2*dz2 +
    double(qq(i    , j + 1, k    ))*dx2*dy1*dz2 +
    double(qq(i    , j    , k + 1))*dx2*dy2*dz1 +
    double(qq(i + 1, j + 1, k    ))*dx1*dy1*dz2 +
    double(qq(i + 1, j    , k + 1))*dx1*dy2*dz1 +
    double(qq(i    , j + 1, k + 1))*dx2*dy1*dz1 +
    double(qq(i + 1, j + 1, k + 1))*dx1*dy1*dz1;
    // clang-format on
}

// Everything the tracer needs, with no Python objects in it, so the loop can
// run with the GIL released and across OpenMP threads.
struct TraceCore
{
    Grid1D gx, gy, gz;
    const ViewArray<float> *bx, *by, *bz;
    std::vector<const ViewArray<float> *> qq;  // interpolated along the line

    bool in_domain(double ox, double oy, double oz) const
    {
        return ox >= 0.0 && ox <= gx.span &&
               oy >= 0.0 && oy <= gy.span &&
               oz >= 0.0 && oz <= gz.span;
    }

    // Unit vector along B.  Returns false where |B| vanishes.
    bool unit_b(double ox, double oy, double oz,
                double &ux, double &uy, double &uz) const
    {
        const double bxt = interpolate3d_offset(*bx, ox, oy, oz, gx, gy, gz);
        const double byt = interpolate3d_offset(*by, ox, oy, oz, gx, gy, gz);
        const double bzt = interpolate3d_offset(*bz, ox, oy, oz, gx, gy, gz);
        const double bb = std::sqrt(bxt * bxt + byt * byt + bzt * bzt);
        if (!(bb > 0.0))
            return false;
        ux = bxt / bb;
        uy = byt / bb;
        uz = bzt / bb;
        return true;
    }

    // Trace one line.  Offsets in, offsets out.  out_x/y/z and each out_q must
    // have room for n_steps entries.  Returns the stop reason and writes the
    // number of populated entries to n_valid.
    int trace_one(double ox0, double oy0, double oz0,
                  double ds, size_t n_steps, double sign,
                  bool implicit, int max_iter, double tol_rel,
                  double *out_x, double *out_y, double *out_z,
                  float *const *out_q,
                  size_t &n_valid) const
    {
        out_x[0] = ox0;
        out_y[0] = oy0;
        out_z[0] = oz0;
        for (size_t f = 0; f < qq.size(); f++)
            out_q[f][0] = (float)interpolate3d_offset(*qq[f], ox0, oy0, oz0, gx, gy, gz);
        n_valid = 1;

        const double h = sign * ds;
        const double tol = tol_rel * ds;

        for (size_t n = 1; n < n_steps; n++)
        {
            const double px = out_x[n - 1], py = out_y[n - 1], pz = out_z[n - 1];
            double ux, uy, uz;
            if (!unit_b(px, py, pz, ux, uy, uz))
                return TRACE_ZERO_FIELD;

            // Predictor, also the answer for the explicit scheme.
            double qx = px + h * ux, qy = py + h * uy, qz = pz + h * uz;

            if (implicit)
            {
                // p_{n+1} = p_n + h * u((p_n + p_{n+1})/2), by fixed point.
                // Time symmetric, so a reversed trace retraces this step.
                bool converged = false;
                for (int it = 0; it < max_iter; it++)
                {
                    const double mx = 0.5 * (px + qx);
                    const double my = 0.5 * (py + qy);
                    const double mz = 0.5 * (pz + qz);
                    if (!in_domain(mx, my, mz))
                        return TRACE_OUT_OF_DOMAIN;
                    if (!unit_b(mx, my, mz, ux, uy, uz))
                        return TRACE_ZERO_FIELD;

                    const double nx = px + h * ux;
                    const double ny = py + h * uy;
                    const double nz = pz + h * uz;
                    const double diff = std::max(std::max(std::abs(nx - qx),
                                                          std::abs(ny - qy)),
                                                 std::abs(nz - qz));
                    qx = nx;
                    qy = ny;
                    qz = nz;
                    if (diff <= tol)
                    {
                        converged = true;
                        break;
                    }
                }
                if (!converged)
                    return TRACE_NO_CONVERGE;
            }
            else
            {
                // Explicit midpoint, the pre-2026-08-23 behaviour.
                const double mx = px + 0.5 * h * ux;
                const double my = py + 0.5 * h * uy;
                const double mz = pz + 0.5 * h * uz;
                if (!in_domain(mx, my, mz))
                    return TRACE_OUT_OF_DOMAIN;
                if (!unit_b(mx, my, mz, ux, uy, uz))
                    return TRACE_ZERO_FIELD;
                qx = px + h * ux;
                qy = py + h * uy;
                qz = pz + h * uz;
            }

            if (!in_domain(qx, qy, qz))
                return TRACE_OUT_OF_DOMAIN;

            out_x[n] = qx;
            out_y[n] = qy;
            out_z[n] = qz;
            for (size_t f = 0; f < qq.size(); f++)
                out_q[f][n] = (float)interpolate3d_offset(*qq[f], qx, qy, qz, gx, gy, gz);
            n_valid = n + 1;
        }
        return TRACE_OK;
    }
};

// Copy the first n_valid elements of a 1D array into a freshly allocated,
// correctly sized array (used to truncate trace output at the point a field
// line left the grid).
inline py::array_t<float> _truncate_1d(const py::array_t<float> &arr_np, size_t n_valid)
{
    py::array_t<float> out_np(n_valid);
    auto in_view = view_array<float>(arr_np);
    auto out_view = view_array<float>(out_np);
    std::copy(in_view.data, in_view.data + n_valid, out_view.data);
    return out_np;
}

inline bool _method_is_implicit(const std::string &method)
{
    if (method == "implicit_midpoint")
        return true;
    if (method == "explicit_midpoint")
        return false;
    throw std::runtime_error(
        "trace_field_line: method must be 'implicit_midpoint' or 'explicit_midpoint'");
}

py::dict trace_field_lines(
    const py::array_t<double, py::array::c_style | py::array::forcecast> &x_np,
    const py::array_t<double, py::array::c_style | py::array::forcecast> &y_np,
    const py::array_t<double, py::array::c_style | py::array::forcecast> &z_np,
    const py::array_t<float, py::array::c_style | py::array::forcecast> &bx_np,
    const py::array_t<float, py::array::c_style | py::array::forcecast> &by_np,
    const py::array_t<float, py::array::c_style | py::array::forcecast> &bz_np,
    py::dict fields,
    const py::array_t<double, py::array::c_style | py::array::forcecast> &x0_np,
    const py::array_t<double, py::array::c_style | py::array::forcecast> &y0_np,
    const py::array_t<double, py::array::c_style | py::array::forcecast> &z0_np,
    double ds, size_t n_steps, double sign = 1.0,
    const std::string &method = "implicit_midpoint",
    int max_iter = 20, double tol = 1e-10)
{
    auto x = view_array<double>(x_np);
    auto y = view_array<double>(y_np);
    auto z = view_array<double>(z_np);

    TraceCore core;
    core.gx = make_grid(x, "x");
    core.gy = make_grid(y, "y");
    core.gz = make_grid(z, "z");

    auto bx = view_array<float>(bx_np);
    auto by = view_array<float>(by_np);
    auto bz = view_array<float>(bz_np);
    core.bx = &bx;
    core.by = &by;
    core.bz = &bz;

    const bool implicit = _method_is_implicit(method);
    if (n_steps < 1)
        throw std::runtime_error("trace_field_lines: n_steps must be >= 1");

    auto s0 = view_array<double>(x0_np);
    auto t0 = view_array<double>(y0_np);
    auto u0 = view_array<double>(z0_np);
    const size_t n_lines = s0.i_size;
    if (t0.i_size != n_lines || u0.i_size != n_lines)
        throw std::runtime_error("trace_field_lines: x0, y0, z0 must have the same length");

    // Seeds, converted to offsets and checked before any tracing starts.
    std::vector<double> ox0(n_lines), oy0(n_lines), oz0(n_lines);
    for (size_t l = 0; l < n_lines; l++)
    {
        ox0[l] = s0(l) - core.gx.origin;
        oy0[l] = t0(l) - core.gy.origin;
        oz0[l] = u0(l) - core.gz.origin;
        if (!core.in_domain(ox0[l], oy0[l], oz0[l]))
        {
            throw std::runtime_error(
                "trace_field_lines: starting point " + std::to_string(l) +
                " is outside the grid domain");
        }
    }

    // Field views, and their output buffers, allocated while we hold the GIL.
    std::vector<std::string> names;
    std::vector<py::array_t<float, py::array::c_style | py::array::forcecast>> qq_np;
    for (auto item : fields)
    {
        names.push_back(py::cast<std::string>(item.first));
        qq_np.push_back(
            py::cast<py::array_t<float, py::array::c_style | py::array::forcecast>>(
                item.second));
    }
    std::vector<ViewArray<float>> qq_views;
    qq_views.reserve(qq_np.size());
    for (auto &a : qq_np)
        qq_views.push_back(view_array<float>(a));
    for (auto &v : qq_views)
        core.qq.push_back(&v);

    const std::vector<py::ssize_t> shape2 = {(py::ssize_t)n_lines, (py::ssize_t)n_steps};
    py::array_t<double> X(shape2), Y(shape2), Z(shape2);
    std::vector<py::array_t<float>> Q;
    Q.reserve(names.size());
    for (size_t f = 0; f < names.size(); f++)
        Q.emplace_back(py::array_t<float>(shape2));

    py::array_t<int64_t> NV((py::ssize_t)n_lines);
    py::array_t<int32_t> ST((py::ssize_t)n_lines);

    double *Xp = X.mutable_data(), *Yp = Y.mutable_data(), *Zp = Z.mutable_data();
    std::vector<float *> Qp;
    for (auto &a : Q)
        Qp.push_back(a.mutable_data());
    int64_t *NVp = NV.mutable_data();
    int32_t *STp = ST.mutable_data();

    {
        py::gil_scoped_release release;
#ifdef _OPENMP
#pragma omp parallel
#endif
        {
            std::vector<float *> q_line(Qp.size());
#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
            for (long long l = 0; l < (long long)n_lines; l++)
            {
                const size_t off = (size_t)l * n_steps;
                for (size_t f = 0; f < Qp.size(); f++)
                    q_line[f] = Qp[f] + off;
                size_t n_valid = 0;
                const int st = core.trace_one(
                    ox0[l], oy0[l], oz0[l], ds, n_steps, sign,
                    implicit, max_iter, tol,
                    Xp + off, Yp + off, Zp + off, q_line.data(), n_valid);
                NVp[l] = (int64_t)n_valid;
                STp[l] = (int32_t)st;
                // Leave the unfilled tail as NaN so a caller cannot silently
                // read stale zeros past the end of a truncated line.
                for (size_t n = n_valid; n < n_steps; n++)
                {
                    Xp[off + n] = std::nan("");
                    Yp[off + n] = std::nan("");
                    Zp[off + n] = std::nan("");
                    for (size_t f = 0; f < Qp.size(); f++)
                        q_line[f][n] = std::nanf("");
                }
            }
        }
    }

    // Back to absolute coordinates.
    for (size_t i = 0; i < n_lines * n_steps; i++)
    {
        Xp[i] += core.gx.origin;
        Yp[i] += core.gy.origin;
        Zp[i] += core.gz.origin;
    }

    py::dict result;
    result["x"] = X;
    result["y"] = Y;
    result["z"] = Z;
    for (size_t f = 0; f < names.size(); f++)
        result[py::str(names[f])] = Q[f];
    result["n_valid"] = NV;
    result["status"] = ST;
    result["ds"] = ds;
    result["method"] = method;
    return result;
}

py::dict trace_field_line(
    const py::array_t<double, py::array::c_style | py::array::forcecast> &x_np,
    const py::array_t<double, py::array::c_style | py::array::forcecast> &y_np,
    const py::array_t<double, py::array::c_style | py::array::forcecast> &z_np,
    const py::array_t<float, py::array::c_style | py::array::forcecast> &bx_np,
    const py::array_t<float, py::array::c_style | py::array::forcecast> &by_np,
    const py::array_t<float, py::array::c_style | py::array::forcecast> &bz_np,
    py::dict fields,
    double x0, double y0, double z0,
    double ds, size_t n_steps, double sign = 1.0,
    const std::string &method = "implicit_midpoint",
    int max_iter = 20, double tol = 1e-10)
{
    auto x = view_array<double>(x_np);
    auto y = view_array<double>(y_np);
    auto z = view_array<double>(z_np);

    TraceCore core;
    core.gx = make_grid(x, "x");
    core.gy = make_grid(y, "y");
    core.gz = make_grid(z, "z");

    auto bx = view_array<float>(bx_np);
    auto by = view_array<float>(by_np);
    auto bz = view_array<float>(bz_np);
    core.bx = &bx;
    core.by = &by;
    core.bz = &bz;

    const bool implicit = _method_is_implicit(method);
    if (n_steps < 1)
        throw std::runtime_error("trace_field_line: n_steps must be >= 1");

    const double ox0 = x0 - core.gx.origin;
    const double oy0 = y0 - core.gy.origin;
    const double oz0 = z0 - core.gz.origin;
    if (!core.in_domain(ox0, oy0, oz0))
        throw std::runtime_error("trace_field_line: starting point is outside the grid domain");

    std::vector<std::string> names;
    std::vector<py::array_t<float, py::array::c_style | py::array::forcecast>> qq_np;
    for (auto item : fields)
    {
        names.push_back(py::cast<std::string>(item.first));
        qq_np.push_back(
            py::cast<py::array_t<float, py::array::c_style | py::array::forcecast>>(
                item.second));
    }
    std::vector<ViewArray<float>> qq_views;
    qq_views.reserve(qq_np.size());
    for (auto &a : qq_np)
        qq_views.push_back(view_array<float>(a));
    for (auto &v : qq_views)
        core.qq.push_back(&v);

    std::vector<double> xb(n_steps), yb(n_steps), zb(n_steps);
    std::vector<std::vector<float>> qb(names.size(), std::vector<float>(n_steps));
    std::vector<float *> qp(names.size());
    for (size_t f = 0; f < names.size(); f++)
        qp[f] = qb[f].data();

    size_t n_valid = 0;
    int st = 0;
    {
        py::gil_scoped_release release;
        st = core.trace_one(ox0, oy0, oz0, ds, n_steps, sign,
                            implicit, max_iter, tol,
                            xb.data(), yb.data(), zb.data(), qp.data(), n_valid);
    }

    py::array_t<double> X((py::ssize_t)n_valid), Y((py::ssize_t)n_valid), Z((py::ssize_t)n_valid);
    double *Xp = X.mutable_data(), *Yp = Y.mutable_data(), *Zp = Z.mutable_data();
    for (size_t n = 0; n < n_valid; n++)
    {
        Xp[n] = xb[n] + core.gx.origin;
        Yp[n] = yb[n] + core.gy.origin;
        Zp[n] = zb[n] + core.gz.origin;
    }

    py::dict result;
    result["x"] = X;
    result["y"] = Y;
    result["z"] = Z;
    for (size_t f = 0; f < names.size(); f++)
    {
        py::array_t<float> a((py::ssize_t)n_valid);
        std::copy(qb[f].data(), qb[f].data() + n_valid, a.mutable_data());
        result[py::str(names[f])] = a;
    }
    result["n_valid"] = (int64_t)n_valid;
    result["status"] = (int32_t)st;
    result["ds"] = ds;
    result["method"] = method;
    return result;
}
