// xtensor-interpolate: https://github.com/rjsberry/xtensor-interpolate
//
// Copyright (C) 2018, Richard Berry <rjsberry@protonmail.com>
//
// Distributed under the terms of BSD 2-Clause "simplified" license. (See
// accompanying file LICENSE, or copy at
// https://github.com/rjsberry/xtensor-interpolate/blob/master/LICENSE)
//

#include "gtest/gtest.h"

#include "xtensor/xmath.hpp"
#include "xtensor/xrandom.hpp"
#include "xtensor/xtensor.hpp"

#include "xtensor-interpolate/xinterpolate.hpp"

namespace xt
{

TEST(univariate_spline, smoke_test)
{
    const std::size_t orig_arr_len = 10;
    const std::size_t test_arr_len = 15;

    std::vector<std::function<xtensor<double, 1> (xtensor<double, 1>)>> tests =
    {
        [](auto x){ return sin(x); },
        [](auto x){ return cos(x); },
        [](auto x){ return x * x * x; },
        [](auto x){ return 1/(1 + 25*x*x); }
    };

    for (const auto& f : tests)
    {
        xtensor<double, 1> x = linspace<double>(-M_PI, M_PI, orig_arr_len);
        xtensor<double, 1> y = f(x);

        auto s = interpolate::UnivariateSpline(x, y).set_smoothing_factor(0.0);

        xtensor<double, 1> x_interp = linspace<double>(-M_PI, M_PI, test_arr_len);
        auto y_interp = s(x_interp);

        EXPECT_TRUE(isclose(y_interp, f(x_interp))());
    }
}

}  // namespace xt
