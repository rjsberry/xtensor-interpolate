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

#include "xtensor-fitpack/xfitpack_utils.hpp"

#include "FCMangle.h"

extern "C" {

void curfit(
        int *iopt, int *m, double *x, double *y, double *w, double *xb,
        double *xe, int *k, double *s, int *nest, int *n, double *t,
        double *c, double *fp, double *wrk, int *lwrk, int *iwrk, int *ier
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
template<class E1, class E2>
auto splrep(xtensor<E1, 1>& x, xt::xtensor<E2, 1>& y, int k = 3) {
    int m = x.size();
    double s = 0.;
    double xb = x[0];
    double xe = x[m-1];
    int iopt = 0;

    // Weights used in computing the weighted least-squares spline fit.
    auto w = allocate_homogeneous_array<double>(m, 1.0);

    int nest = std::max(m + k + 1, 2*k + 3);

    // Knots.
    auto t = allocate_homogeneous_array<double>(nest, 0.0);

    // Coefficients.
    auto c = allocate_homogeneous_array<double>(nest, 0.0);

    // Working memory.
    int lwrk = m*(k + 1) + nest*(7 + 3*k);
    auto wrk = allocate_homogeneous_array<double>(lwrk, 0.0);
    auto iwrk = allocate_homogeneous_array<int>(nest, 0);

    double fp = 0.0;
    int n = 0;
    int ier = 0;

    curfit(
        &iopt, &m, &x[0], &y[0], w.get(), &xb, &xe, &k, &s, &nest, &n,
        t.get(), c.get(), &fp, wrk.get(), &lwrk, iwrk.get(), &ier
    );

    trim_array_ptr(t, nest, n);
    trim_array_ptr(c, nest, n);

    std::vector<std::size_t> shape = { static_cast<std::size_t>(n) };
    xtensor<double, 1> t_out = adapt(t.get(), n, no_ownership(), shape);
    xtensor<double, 1> c_out = adapt(t.get(), n, no_ownership(), shape);

    auto tck = std::make_tuple(t_out, c_out, k);

    return tck;
}

}  // fitpack
}  // namespace xt

#endif  // INCLUDE_XTENSOR_FITPACK_XFITPACK_HPP_
