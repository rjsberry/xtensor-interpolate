#ifndef INCLUDE_XTENSOR_INTERPOLATE_XINTERPOLATE_HPP_
#define INCLUDE_XTENSOR_INTERPOLATE_XINTERPOLATE_HPP_

/// xtensor-interpolate: https://github.com/rjsberry/xtensor-interpolate
///
/// Copyright (C) 2018, Richard Berry <rjsberry@protonmail.com>
///
/// Distributed under the terms of BSD 2-Clause "simplified" license. (See
/// accompanying file LICENSE, or copy at
/// https://github.com/rjsberry/xtensor-interpolate/blob/master/LICENSE)
///

#include <map>
#include <string>
#include <tuple>
#include <vector>

#include "xtensor/xarray.hpp"
#include "xtensor/xexpression.hpp"
#include "xtensor/xlayout.hpp"
#include "xtensor/xtensor.hpp"
#include "xtensor/xview.hpp"

#include "xfitpack/xfcdec.h"

namespace xt
{
namespace interpolate
{

/// Find the B-spline representation of a 1-D curve.
///
/// Given the set of data points `(x[i], y[i])` determine a smooth spline
/// approximation of degree k on the interval `xb <= x <= xe`.
///
/// @param [in] x,y
///     The data points defining a curve y = f(x).
/// @param [in] k
///     The order of the spline fit. It is recommended to use cubic splines.
///     1 <= k <= 5
///
/// @returns A tuple containing the vector of knots, `t`, the B-spline
///          coefficients, `c`, and the degree of the spline, `k`. It is not
///          recommended to manually adjust any of these values.
///
/// @throws std::runtime_error
///     Thrown if ``FITPACK`` encounters any errors.
///
/// @todo Address spline classification based on ``ier`` return.
///
/// @todo Catch and solve ``nest`` size error.
///
template <class E1, class E2>
auto splrep(const xexpression<E1>& x, const xexpression<E2>& y, int k = 3)
{
    auto m = static_cast<int>(x.derived_cast().shape()[0]);

    auto s = 0.0;
    auto xb = x.derived_cast()[0];
    auto xe = x.derived_cast()[m - 1];
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

    _fc_curfit(&iopt, &m, &x.derived_cast()[0], &y.derived_cast()[0], &w[0],
               &xb, &xe, &k, &s, &nest, &n, &t[0], &c[0], &fp, &wrk[0], &lwrk,
               &iwrk[0], &ier);

    switch (ier)
    {
      case -2:
      case -1:
      case 0:
        break;  // Normal return.
      case 1:
        throw std::runtime_error("recalculation of nest is required");
      case 2:
        throw std::runtime_error("theoritcally impossible result");
      case 3:
        throw std::runtime_error("reached maximum iterations");
      case 10:
        throw std::runtime_error("invalid input data");
      default:
        throw std::runtime_error("an unknown error occurred");
    }

    xtensor<double, 1> _t = view(t, range(0, n));
    xtensor<double, 1> _c = view(c, range(0, n));

    return std::make_tuple(_t, _c, k);
}

/// Evaluate a B-Spline or its derivatives.
///
/// @param [in] x
///     An array of points at which to return the value of the smoothed
///     spline or its derivatives.
/// @param [in] tck
///     The tuple returned by `splrep` containing the knots, coefficients,
///     and degree of the spline.
/// @param [in] der
///     The order of derivative of the spline to compute (must be less than
///      or equal to k).
/// @param [in] ext
///     Controls the value returned for elements of `x` not in the interval
///     defined by the knot sequence.
///       * if `0`, return the extrapolated value.
///       * if `1`, return 0
///       * if `2`, throw exception.
///       * if `3`, return the boundary value.
///
/// @returns An array of interpolated values corresponding to the input
///          parameter `x`.
///
/// @throws std::runtime_error
///     Thrown if ``FITPACK`` encounters any errors.
///
template <class E, class... Args>
auto splev(const xexpression<E>& x,
           const std::tuple<Args...>& tck,
           int der = 0,
           int ext = 0)
{
    auto m = static_cast<int>(x.derived_cast().shape()[0]);
    auto n = static_cast<int>(std::get<0>(tck).size());

    xtensor<double, 1> y = zeros<double>({ static_cast<std::size_t>(m) });
    auto ier = 0;

    _fc_splev(&std::get<0>(tck)[0], &n, &std::get<1>(tck)[0], &std::get<2>(tck),
              &x.derived_cast()[0], &y[0], &m, &ext, &ier);

    switch (ier)
    {
      case 0:
        break;  // Normal return
      case 1:
        throw std::runtime_error("evaluation outside bounds of the spline");
      case 10:
        throw std::runtime_error("invalid input data");
      default:
        throw std::runtime_error("an unknown error occurred");
    }

    return y;
}

/// Evaluate the definite integral of a B-spline between two points.
///
/// @param [in] a,b
///     The endpoints defining the bounds of integration.
/// @param [in] tck
///     A tuple containing the knots, B-spline coefficients, and degree of
///     the spline.
///
/// @returns The resultant integral.
///
/// @note This routine assums the spline is **0** outside of its data points.
///
template <class... Args>
auto splint(double a, double b, const std::tuple<Args...>& tck)
{
    auto n = static_cast<int>(std::get<0>(tck).size());

    xtensor<double, 1> wrk = zeros<double>({ n });

    return _fc_splint(&std::get<0>(tck)[0], &n, &std::get<1>(tck)[0],
                      &std::get<2>(tck), &a, &b, &wrk[0]);
}

/// Evaluate all derivatives of a B-spline.
///
/// @param [in] x
///     A set of points at which to evaluate the derivatives.
/// @param [in] tck
///     A tuple containing the knots, B-spline coefficients, and degree of
///     the spline.
///
/// @returns An `xarray` of all derivatives of the spline up to order `k` for
///          each point in `x`.
///
/// @throws std::runtime_error
///     Thrown if ``FITPACK`` encounters any errors.
///
template <class E, class... Args>
auto spalde(const xexpression<E>& x, const std::tuple<Args...>& tck)
{
    auto m = x.derived_cast().shape()[0];
    auto n = static_cast<int>(std::get<0>(tck).size());

    auto k1 = std::get<2>(tck) + 1;

    xarray<double> d = zeros<double>({ static_cast<std::size_t>(k1) * m });
    for (std::size_t i = 0; i < m; ++i)
    {
        auto ier = 0;
        _fc_spalde(&std::get<0>(tck)[0], &n, &std::get<1>(tck)[0], &k1,
                   &x.derived_cast()[i], &d[i * k1], &ier);
                   
        switch (ier)
        {
        case 0:
            break;  // Normal return
        case 10:
            throw std::runtime_error("invalid input data");
        default:
            throw std::runtime_error("an unknown error occurred");
        }
    }
    d.reshape({ m, static_cast<std::size_t>(k1) });

    return d;
}

/// Compute the spline representation of the derivative of a given spline.
///
/// @param [in] x
///     An array of points at which to evaluate the derivative.
/// @param [in] tck
///     A tuple containing the knots, B-spline coefficients, and degree of
///     the spline.
/// @param [in] nu
///     Order of derivative to evaluate.
///
/// @returns A tuple containing the vector of knots, `t`, the B-spline
///          coefficients, `c`, and the degree of the spline, `k - nu` for the
///          evaluated derivative. It is not recommended to manually adjust any 
///          of these values.
///
/// @throws std::invalid_argument
///     Thrown if `0 <= n <= k` does not hold, where `k` is the order of the
///     spline.
/// @throws std::runtime_error
///     Thrown if ``FITPACK`` encounters any errors.
///
/// @todo Link to ``splantider`` when `n < 0`.
///
/// @todo Parametrize ``ext`` to allow users to control FITPACK behaviour.
///
template <class E, class... Args>
auto splder(const xexpression<E>& x, const std::tuple<Args...>& tck, int nu = 1)
{
    if (nu > std::get<2>(tck))
    {
        throw std::invalid_argument("order of derivative must be <= k");
    }
    if (nu < 0)
    {
        throw std::invalid_argument("not implemented");
    }

    auto n = static_cast<int>(std::get<0>(tck).shape()[0]);
    auto m = static_cast<int>(x.derived_cast().shape()[0]);

    xtensor<double, 1> y = zeros<double>({ x.derived_cast().shape()[0] });
    xtensor<double, 1> wrk = zeros<double>({ n });
    auto ext = 0;
    auto ier = 0;

    _fc_splder(&std::get<0>(tck)[0], &n, &std::get<1>(tck)[0], &std::get<2>(tck),
               &nu, &x.derived_cast()[0], &y[0], &m, &ext, &wrk[0], &ier);

    switch (ier)
    {
      case 0:
        break;  // Normal return
      case 1:
        throw std::runtime_error("evaluation outside bounds of the spline");
      case 10:
        throw std::runtime_error("invalid input data");
      default:
        throw std::runtime_error("an unknown error occurred");
    }

    return splrep(x, y, std::get<2>(tck) - nu);
}

}  // namespace interpolate

}  // namespace xt

#endif  // INCLUDE_XTENSOR_INTERPOLATE_XINTERPOLATE_HPP_
