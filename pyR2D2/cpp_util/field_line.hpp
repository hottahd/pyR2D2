#pragma once
#include <cmath>
#include <omp.h>
#include "view_array.hpp"

// This file contains the implementation of field line tracing and interpolation.
// The main function is trace_field_line
struct FieldData
{
    std::string name;

    py::array_t<float, py::array::c_style | py::array::forcecast> qq_np;
    py::array_t<float> out_np;

    ViewArray<float> qq;
    ViewArray<float> out;
};

inline float interpolate3d(const ViewArray<float> &qq, float xt, float yt, float zt, const ViewArray<float> &x, const ViewArray<float> &y, const ViewArray<float> &z)
{
    float dx = x(1) - x(0);
    float dy = y(1) - y(0);
    float dz = z(1) - z(0);
    int i = (int)std::floor((xt - x(0)) / dx);
    int j = (int)std::floor((yt - y(0)) / dy);
    int k = (int)std::floor((zt - z(0)) / dz);

    float dx1 = (xt - x(i)) / dx;
    float dy1 = (yt - y(j)) / dy;
    float dz1 = (zt - z(k)) / dz;
    float dx2 = 1.0 - dx1;
    float dy2 = 1.0 - dy1;
    float dz2 = 1.0 - dz1;

    // clang-format off
    return 
    qq(i    , j    , k    )*dx2*dy2*dz2 +
    qq(i + 1, j    , k    )*dx1*dy2*dz2 + 
    qq(i    , j + 1, k    )*dx2*dy1*dz2 +
    qq(i    , j    , k + 1)*dx2*dy2*dz1 +
    qq(i + 1, j + 1, k    )*dx1*dy1*dz2 + 
    qq(i + 1, j    , k + 1)*dx1*dy2*dz1 +
    qq(i    , j + 1, k + 1)*dx2*dy1*dz1 +
    qq(i + 1, j + 1, k + 1)*dx1*dy1*dz1;
    // clang-format on
}

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

py::dict trace_field_line(
    const py::array_t<float, py::array::c_style | py::array::forcecast> &x_np,
    const py::array_t<float, py::array::c_style | py::array::forcecast> &y_np,
    const py::array_t<float, py::array::c_style | py::array::forcecast> &z_np,
    const py::array_t<float, py::array::c_style | py::array::forcecast> &bx_np,
    const py::array_t<float, py::array::c_style | py::array::forcecast> &by_np,
    const py::array_t<float, py::array::c_style | py::array::forcecast> &bz_np,
    py::dict fields,
    float x0, float y0, float z0,
    float ds, size_t n_steps, float sign = 1.0)
{
    auto x = view_array<float>(x_np);
    auto y = view_array<float>(y_np);
    auto z = view_array<float>(z_np);

    // Grids are assumed monotonically increasing, matching the assumption
    // already made elsewhere in this codebase (e.g. eos.hpp, rte.hpp).
    // The upper bound is tightened by one cell so that interpolate3d's
    // (i, j, k) and (i+1, j+1, k+1) accesses always stay in range.
    float xmin = x(0), xmax = x(x.i_size - 2);
    float ymin = y(0), ymax = y(y.i_size - 2);
    float zmin = z(0), zmax = z(z.i_size - 2);

    auto in_domain = [&](float xt, float yt, float zt) {
        return xt >= xmin && xt <= xmax &&
               yt >= ymin && yt <= ymax &&
               zt >= zmin && zt <= zmax;
    };

    if (!in_domain(x0, y0, z0))
    {
        throw std::runtime_error("trace_field_line: starting point is outside the grid domain");
    }

    py::array_t<float> x_fl_np(n_steps), y_fl_np(n_steps), z_fl_np(n_steps);
    auto x_fl = view_array<float>(x_fl_np);
    auto y_fl = view_array<float>(y_fl_np);
    auto z_fl = view_array<float>(z_fl_np);

    x_fl(0) = x0;
    y_fl(0) = y0;
    z_fl(0) = z0;

    std::vector<FieldData> field_datas;
    field_datas.reserve(fields.size());

    for (auto item : fields)
    {
        FieldData fd;
        fd.name = py::cast<std::string>(item.first);
        fd.qq_np = py::cast<py::array_t<float, py::array::c_style | py::array::forcecast>>(item.second);
        fd.out_np = py::array_t<float>(n_steps);
        fd.qq = view_array<float>(fd.qq_np);
        fd.out = view_array<float>(fd.out_np);

        field_datas.emplace_back(std::move(fd));
    }

    auto bx = view_array<float>(bx_np);
    auto by = view_array<float>(by_np);
    auto bz = view_array<float>(bz_np);

    for (auto &fd : field_datas)
    {
        fd.out(0) = interpolate3d(fd.qq, x_fl(0), y_fl(0), z_fl(0), x, y, z);
    }

    // Number of leading entries in x_fl/y_fl/z_fl/fd.out that are actually
    // populated. Tracing stops as soon as the field line (or the half-step
    // point used to advance it) leaves the grid domain, since B is not
    // defined -- and therefore not safe to read -- outside it.
    size_t n_valid = 1;
    {
        py::gil_scoped_release release;
        for (size_t n = 1; n < n_steps; n++)
        {
            float bxt = interpolate3d(bx, x_fl(n - 1), y_fl(n - 1), z_fl(n - 1), x, y, z);
            float byt = interpolate3d(by, x_fl(n - 1), y_fl(n - 1), z_fl(n - 1), x, y, z);
            float bzt = interpolate3d(bz, x_fl(n - 1), y_fl(n - 1), z_fl(n - 1), x, y, z);
            float bb = std::sqrt(bxt * bxt + byt * byt + bzt * bzt);

            float xm_fl = x_fl(n - 1) + 0.5 * sign * ds * bxt / bb;
            float ym_fl = y_fl(n - 1) + 0.5 * sign * ds * byt / bb;
            float zm_fl = z_fl(n - 1) + 0.5 * sign * ds * bzt / bb;

            if (!in_domain(xm_fl, ym_fl, zm_fl))
                break;

            bxt = interpolate3d(bx, xm_fl, ym_fl, zm_fl, x, y, z);
            byt = interpolate3d(by, xm_fl, ym_fl, zm_fl, x, y, z);
            bzt = interpolate3d(bz, xm_fl, ym_fl, zm_fl, x, y, z);
            bb = std::sqrt(bxt * bxt + byt * byt + bzt * bzt);

            float x_new = x_fl(n - 1) + sign * ds * bxt / bb;
            float y_new = y_fl(n - 1) + sign * ds * byt / bb;
            float z_new = z_fl(n - 1) + sign * ds * bzt / bb;

            if (!in_domain(x_new, y_new, z_new))
                break;

            x_fl(n) = x_new;
            y_fl(n) = y_new;
            z_fl(n) = z_new;

            for (auto &fd : field_datas)
            {
                fd.out(n) = interpolate3d(fd.qq, x_fl(n), y_fl(n), z_fl(n), x, y, z);
            }

            n_valid = n + 1;
        }
    }

    py::dict result;
    result["x"] = _truncate_1d(x_fl_np, n_valid);
    result["y"] = _truncate_1d(y_fl_np, n_valid);
    result["z"] = _truncate_1d(z_fl_np, n_valid);
    for (auto &fd : field_datas)
    {
        result[py::str(fd.name)] = _truncate_1d(fd.out_np, n_valid);
    }
    return result;
}