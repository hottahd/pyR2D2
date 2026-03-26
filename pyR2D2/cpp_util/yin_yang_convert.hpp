#pragma once
#include <array>
#include <iostream>
#include <cmath>
#include <omp.h>
#include "const.hpp"
#include "view_array.hpp"

// Lagrange interpolation for 3rd order (4 points)
// x: array of 4 points, xt: target point
template <typename T>
inline std::array<T, 4> lagrange_interpolation_weight_3rd(const std::array<T, 4> &x, T xt)
{
    std::array<T, 4> weight;
    weight[0] = ((xt - x[1]) * (xt - x[2]) * (xt - x[3])) / ((x[0] - x[1]) * (x[0] - x[2]) * (x[0] - x[3]));
    weight[1] = ((xt - x[0]) * (xt - x[2]) * (xt - x[3])) / ((x[1] - x[0]) * (x[1] - x[2]) * (x[1] - x[3]));
    weight[2] = ((xt - x[0]) * (xt - x[1]) * (xt - x[3])) / ((x[2] - x[0]) * (x[2] - x[1]) * (x[2] - x[3]));
    weight[3] = ((xt - x[0]) * (xt - x[1]) * (xt - x[2])) / ((x[3] - x[0]) * (x[3] - x[1]) * (x[3] - x[2]));

    return weight;
}

// 2D Lagrange interpolation for 3rd order (4 points in each direction)
// tht: target theta, dth_yy_o: inverse of theta grid spacing in Yin-Yang grid
// pht: target phi, dph_yy_o: inverse of phi grid spacing in
// Yin-Yang grid, th_yy: theta grid points in Yin-Yang grid, ph_yy: phi grid points in Yin-Yang grid
// qq_yy: 3D array of the quantity to be interpolated in Yin-Yang grid (i, j, k)
// Returns the interpolated value at (tht, pht) in regular spherical geometry
template <typename T>
inline std::vector<T> lagrange_interpolation_3rd(double tht, double dth_yy_o, double pht, double dph_yy_o, auto const &th_yy, auto const &ph_yy, auto const &qq_yy, const size_t k_size)
{
    size_t ic = (tht - th_yy(0)) * dth_yy_o;
    size_t jc = (pht - ph_yy(0)) * dph_yy_o;

    std::array<size_t, 4> i_s = {static_cast<size_t>(ic - 1), static_cast<size_t>(ic), static_cast<size_t>(ic + 1), static_cast<size_t>(ic + 2)};
    std::array<double, 4> th_vec = {th_yy(i_s[0]), th_yy(i_s[1]), th_yy(i_s[2]), th_yy(i_s[3])};

    std::array<size_t, 4> j_s = {static_cast<size_t>(jc - 1), static_cast<size_t>(jc), static_cast<size_t>(jc + 1), static_cast<size_t>(jc + 2)};
    std::array<double, 4> ph_vec = {ph_yy(j_s[0]), ph_yy(j_s[1]), ph_yy(j_s[2]), ph_yy(j_s[3])};

    auto weight_th = lagrange_interpolation_weight_3rd<double>(th_vec, tht);
    auto weight_ph = lagrange_interpolation_weight_3rd<double>(ph_vec, pht);

    std::vector<T> qq_output_core(k_size);
    for (size_t k = 0; k < k_size; ++k)
    {
        qq_output_core[k] = T(0);
    }

    for (int m = 0; m < 4; ++m)
    {
        for (int n = 0; n < 4; ++n)
        {
            for (size_t k = 0; k < k_size; ++k)
            {
                qq_output_core[k] += weight_th[m] * weight_ph[n] * qq_yy(i_s[m], j_s[n], k);
            }
        }
    }
    return qq_output_core;
}

struct YinYang
{
    py::array_t<double, py::array::c_style | py::array::forcecast> th_yy_np;
    py::array_t<double, py::array::c_style | py::array::forcecast> ph_yy_np;
    py::array_t<double, py::array::c_style | py::array::forcecast> th_np;
    py::array_t<double, py::array::c_style | py::array::forcecast> ph_np;

    YinYang(
        py::array_t<double, py::array::c_style | py::array::forcecast> th_yy_np_,
        py::array_t<double, py::array::c_style | py::array::forcecast> ph_yy_np_,
        py::array_t<double, py::array::c_style | py::array::forcecast> th_np_,
        py::array_t<double, py::array::c_style | py::array::forcecast> ph_np_)
        : th_yy_np(th_yy_np_), ph_yy_np(ph_yy_np_), th_np(th_np_), ph_np(ph_np_) {
          };

    template <typename T>
    py::object convert_scalar(
        const py::array_t<T, py::array::c_style | py::array::forcecast> &qq_yin_np,
        const py::array_t<T, py::array::c_style | py::array::forcecast> &qq_yan_np)
    {

        auto qq_yin = view_array<T>(qq_yin_np);
        auto qq_yan = view_array<T>(qq_yan_np);
        auto th_yy = view_array<double>(th_yy_np);
        auto ph_yy = view_array<double>(ph_yy_np);

        auto th = view_array<double>(th_np);
        auto ph = view_array<double>(ph_np);

        size_t margin = 2;
        size_t i_size = (qq_yin.i_size - 2 * margin) * 2;
        size_t j_size = (qq_yin.i_size - 2 * margin) * 4;
        size_t k_size = qq_yin.k_size;

        double th_min = th_yy(1);
        double th_max = th_yy(th_yy.size - 3);

        double ph_min = ph_yy(1);
        double ph_max = ph_yy(ph_yy.size - 3);

        double dth_yy_o = 1.0 / (th_yy(1) - th_yy(0));
        double dph_yy_o = 1.0 / (ph_yy(1) - ph_yy(0));

        // Output array in regular spherical geometry
        py::array_t<T, py::array::c_style> qq_output_np({i_size, j_size, k_size});
        auto qq_output = view_array<T>(qq_output_np);

        std::vector<T> qq_output_core(k_size);
        // for Yang grid interpolation
        for (size_t i = 0; i < qq_output.i_size; ++i)
        {
            for (size_t j = 0; j < qq_output.j_size; ++j)
            {
                // for Yin grid
                if (th_min <= th(i) && th(i) <= th_max && ph_min <= ph(j) && ph(j) <= ph_max)
                {
                    qq_output_core = lagrange_interpolation_3rd<T>(th(i), dth_yy_o, ph(j), dph_yy_o, th_yy, ph_yy, qq_yin, k_size);
                    for (size_t k = 0; k < k_size; ++k)
                    {
                        qq_output(i, j, k) = qq_output_core[k];
                    }
                }
                else
                {
                    double tht = std::acos(std::sin(th(i)) * std::sin(ph(j)));
                    double pht = std::atan2(std::cos(th(i)), -std::sin(th(i)) * std::cos(ph(j)));

                    qq_output_core = lagrange_interpolation_3rd<T>(tht, dth_yy_o, pht, dph_yy_o, th_yy, ph_yy, qq_yan, k_size);
                    for (size_t k = 0; k < k_size; ++k)
                    {
                        qq_output(i, j, k) = qq_output_core[k];
                    }
                }
            }
        }

        if (qq_yin.ndim == 2)
        {
            return qq_output_np.attr("reshape")(py::make_tuple(i_size, j_size));
        }
        else if (qq_yin.ndim == 3)
        {
            return qq_output_np;
        }
    }
};
