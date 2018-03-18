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

#include "xtensor/xtensor.hpp"
#include "xtensor/xview.hpp"

#include "xfitpack/xfcdec.h"

namespace xt
{
namespace interpolate
{

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
template <class E, class... Args>
auto splev(const xexpression<E>& x,
           const std::tuple<Args...>& tck,
           int der = 0,
           int ext = 0)
{
    auto _x = x.derived_cast();
    auto m = static_cast<int>(x.derived_cast().shape()[0]);

    auto t = std::get<0>(tck);
    auto n = static_cast<int>(t.size());
    auto c = std::get<1>(tck);
    auto k = std::get<2>(tck);

    xtensor<double, 1> y = zeros<double>({ static_cast<std::size_t>(m) });
    auto ier = 0;

    _fc_splev(&t[0], &n, &c[0], &k, &_x[0], &y[0], &m, &ext, &ier);

    return y;
}

}  // namespace interpolate

}  // namespace xt

#endif  // INCLUDE_XTENSOR_INTERPOLATE_XINTERPOLATE_HPP_
