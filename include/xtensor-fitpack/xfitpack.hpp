#ifndef INCLUDE_XTENSOR_FITPACK_XFITPACK_HPP_
#define INCLUDE_XTENSOR_FITPACK_XFITPACK_HPP_

// xtensor-fitpack: https://github.com/rjsberry/xtensor-fitpack
//
// Copyright (C) 2018, Richard Berry <rjsberry@protonmail.com>
//
// Distributed under the terms of BSD 2-Clause "simplified" license. (See
// accompanying file LICENSE, or copy at
// https://github.com/rjsberry/xtensor-fitpack/blob/master/LICENSE)
//

#include <tuple>

#include "xtensor/xexpression.hpp"
#include "xtensor/xtensor.hpp"
#include "xtensor/xview.hpp"

#include "xfcdec.h"

namespace xt
{
namespace interpolate
{

// Find the B-spline representation of a 1-D curve.
//
// Given the set of data points `(x[i], y[i]) determine a smooth spline
// approximation of degree k on the interval `xb <= x <= xe`.
//
// @param [in] x, y
//     The data points defining a curve y = f(x).
// @param [in] k
//     The order of the spline fit. It is recommended to use cubic splines.
//     1 <= k <= 5
//
template <class E1, class E2>
auto splrep(xexpression<E1>& x, xexpression<E2>& y, int k = 3)
{
    XTENSOR_ASSERT(x.dimension() == 1);
    XTENSOR_ASSERT(y.dimension() == 1);

    auto _x = x.derived_cast();
    auto _y = y.derived_cast();
    auto m = static_cast<int>(_x.shape()[0]);

    auto s = 0.0;
    auto xb = _x[0];
    auto xe = _x[m-1];
    auto iopt = 0;

    // Weights used in computing the weighted least-squares spline fit.
    xtensor<double, 1> w = ones<double>({ static_cast<std::size_t>(m) });

    auto nest = std::max(m + k + 1, 2*k + 3);
    // Knots.
    xtensor<double, 1> t = zeros<double>({ static_cast<std::size_t>(nest) });
    // Coefficients.
    xtensor<double, 1> c = zeros<double>({ static_cast<std::size_t>(nest) });

    // Working memory.
    auto lwrk = m*(k + 1) + nest*(7 + 3*k);
    xtensor<double, 1> wrk = zeros<double>({ static_cast<std::size_t>(lwrk) });
    xtensor<int, 1> iwrk = zeros<int>({ static_cast<std::size_t>(lwrk) });

    auto fp = 0.0;
    auto n = 0;
    auto ier = 0;

    _fc_curfit(&iopt, &m, &_x[0], &_y[0], &w[0], &xb, &xe, &k, &s, &nest, &n,
               &t[0], &c[0], &fp, &wrk[0], &lwrk, &iwrk[0], &ier);

    xtensor<double, 1> _t = view(t, range(0, n));
    xtensor<double, 1> _c = view(c, range(0, n));

    return std::make_tuple(_t, _c, k);
}

// Evaluate a B-Spline or its derivatives.
//
// @param [in] x
//     An array of points at which to return the value of the smoothed
//     spline or its derivatives.
// @param [in] tck
//     The tuple returned by `splrep` containing the knots, coefficients,
//     and degree of the spline.
// @param [in] der
//     The order of derivative of the spline to compute (must be less than
//      or equal to k).
// @param [in] ext
//     Controls the value returned for elements of `x` not in the interval
//     defined by the knot sequence.
//       * if ext=0, return the extrapolated value.
//       * if ext=1, return 0
//       * if ext=2, throw `out_of_range`
//       * if ext=3, return the boundary value.
//
template <class E, class... Args>
auto splev(xexpression<E>& x, std::tuple<Args...>& tck, int der = 0, int ext = 0)
{
    auto _x = x.derived_cast();
    auto m = static_cast<int>(x.derived_cast().shape()[0]);

    auto t = std::get<0>(tck);
    auto n = static_cast<int>(t.size());
    auto c = std::get<1>(tck);
    auto k = std::get<2>(tck);

    if (!(0 <= der && der <= k))
    {
        throw std::out_of_range("der");
    }
    if (!(0 <= ext && ext <= 3))
    {
        throw std::out_of_range("ext");
    }

    xtensor<double, 1> y = zeros<double>({ static_cast<std::size_t>(m) });
    auto ier = 0;

    _fc_splev(&t[0], &n, &c[0], &k, &_x[0], &y[0], &m, &ext, &ier);

    return y;
}

}  // namespace interpolate

}  // namespace xt

#endif  // INCLUDE_XTENSOR_FITPACK_XFITPACK_HPP_
