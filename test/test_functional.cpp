// xtensor-interpolate: https://github.com/rjsberry/xtensor-interpolate
//
// Copyright (C) 2018, Richard Berry <rjsberry@protonmail.com>
//
// Distributed under the terms of BSD 2-Clause "simplified" license. (See
// accompanying file LICENSE, or copy at
// https://github.com/rjsberry/xtensor-interpolate/blob/master/LICENSE)
//

#include "gtest/gtest.h"

#define XTENSOR_ENABLE_ASSERT
#include "xtensor/xexception.hpp"
#include "xtensor/xmath.hpp"
#include "xtensor/xtensor.hpp"

#include "xtensor-interpolate/xinterpolate.hpp"

namespace xt
{

TEST(splrep_splev, smoke_test)
{
    xtensor<double, 1> x = linspace<double>(0, M_PI, 10);
    xtensor<double, 1> y = sin(x);

    auto tck = interpolate::splrep(x, y);

    xtensor<double, 1> xs = arange<double>(0, M_PI, 100);
    xtensor<double, 1> ys = interpolate::splev(xs, tck);

    EXPECT_TRUE(allclose(ys, sin(xs), pow10(-5)));
}

}  // namespace xt
