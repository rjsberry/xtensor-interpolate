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

#include <algorithm>
#include <memory>
#include <tuple>

#include "xtensor/xadapt.hpp"
#include "xtensor/xexpression.hpp"
#include "xtensor/xtensor.hpp"

#include "FCMangle.h"

extern "C" {

void _curfit(
        int *iopt, int *m, double *x, double *y, double *w, double *xb,
        double *xe, int *k, double *s, int *nest, int *n, double *t,
        double *c, double *fp, double *wrk, int *lwrk, int *iwrk, int *ier
    );

void _splev(
        double *t, int *n, double *c, int *k, double *x, double *y, int *m,
        int *e, int *ier
    );

}  // extern "C"

namespace xt {
namespace fitpack {

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
auto splrep(xexpression<E1>& x, xexpression<E2>& y, int k = 3) {
    auto m = static_cast<int>(x.derived_cast().size());

    auto s = 0.0;
    auto xb = x.derived_cast()[0];
    auto xe = x.derived_cast()[m-1];
    auto iopt = 0;

    // Weights used in computing the weighted least-squares spline fit.
    std::vector<double> w(static_cast<std::size_t>(m), 1.0);

    auto nest = std::max(m + k + 1, 2*k + 3);

    // Knots.
    std::vector<double> t(static_cast<std::size_t>(nest), 0.0);
    // Coefficients.
    std::vector<double> c(static_cast<std::size_t>(nest), 0.0);
    // Working memory.
    auto lwrk = m*(k + 1) + nest*(7 + 3*k);
    std::vector<double> wrk(static_cast<std::size_t>(lwrk), 0.0);
    std::vector<int> iwrk(static_cast<std::size_t>(nest), 0.0);

    auto fp = 0.0;
    auto n = 0;
    auto ier = 0;

    _curfit(
        &iopt, &m, &x.derived_cast()[0], &y.derived_cast()[0], &w[0], &xb, &xe,
        &k, &s, &nest, &n, &t[0], &c[0], &fp, &wrk[0], &lwrk, &iwrk[0], &ier
    );

    std::vector<std::size_t> tc_shape = { static_cast<std::size_t>(n) };
    auto tx = xt::adapt(t, tc_shape);
    auto cx = xt::adapt(c, tc_shape);
    auto tck = std::make_tuple(tx, cx, k);

    return tck;
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
template <class E1, class E2, class E3>
auto splev(xexpression<E1>& x,
           std::tuple<xexpression<E2>, xexpression<E3>, int>& tck,
           int der = 0,
           int ext = 0) {
    auto m = static_cast<int>(x.size());

    auto t = std::get<0>(tck);
    auto n = static_cast<int>(t.size());
    auto c = std::get<1>(tck);
    auto k = std::get<2>(tck);

    if (!(0 <= der && der <= k)) {
        throw std::out_of_range("der");
    }
    if (!(0 <= ext && ext <= 3)) {
        throw std::out_of_range("ext");
    }

    std::vector<double> y(static_cast<std::size_t>(m), 0);
    int ier = 0;
    _splev(
        &t[0], &n, &c[0], &k, &x[0], &y[0], &m, &ext, &ier
    );

    return y;
}

}  // fitpack
}  // namespace xt

#endif  // INCLUDE_XTENSOR_FITPACK_XFITPACK_HPP_
