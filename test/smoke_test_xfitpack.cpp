// xtensor-fitpack: https://github.com/rjsberry/xtensor-fitpack
//
// Copyright (C) 2018, Richard Berry <rjsberry@protonmail.com>
//
// Distributed under the terms of BSD 2-Clause "simplified" license. (See
// accompanying file LICENSE, or copy at
// https://github.com/rjsberry/xtensor-fitpack/blob/master/LICENSE)
//

// Smoke tests for fitpack routines: generally just check they run safely.

#include <cmath>
#include <vector>

#include "gtest/gtest.h"

#include "xtensor/xmath.hpp"
#include "xtensor/xtensor.hpp"
#include "xtensor/xrandom.hpp"

#include "xtensor-fitpack/xfitpack.hpp"

namespace xt
{

TEST(smoke_test, cubic_bspline_interpolation)
{
    const std::size_t orig_arr_len = 100;
    const std::size_t test_arr_len = 1000;

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

        auto tck = interpolate::splrep(x, y);

        xtensor<double, 1> x_interp = linspace<double>(-M_PI, M_PI, test_arr_len);
        auto y_interp = interpolate::splev(x_interp, tck);

        EXPECT_TRUE(isclose(y_interp, f(x_interp))());
    }
}

}  // namespace xt
