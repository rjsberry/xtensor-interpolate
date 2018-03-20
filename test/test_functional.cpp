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
#include "xtensor/xstrided_view.hpp"
#include "xtensor/xtensor.hpp"

#include "xtensor-interpolate/xinterpolate.hpp"

namespace xt
{

namespace
{

const std::size_t ARR_LEN = 100;

}  // namespace

TEST(smoke_test, splrep_splev)
{
    xtensor<double, 1> x = linspace<double>(0, M_PI, 10);
    xtensor<double, 1> y = sin(x);

    auto tck = interpolate::splrep(x, y);

    xtensor<double, 1> xs = arange<double>(0, M_PI, 100);
    xtensor<double, 1> ys = interpolate::splev(xs, tck);

    EXPECT_TRUE(allclose(ys, sin(xs), pow10(-5)));
}

TEST(smoke_test, splrep_splint)
{
    struct TestCase
    {
        xtensor<double, 1>                                    x;
        std::function<xtensor<double, 1>(xtensor<double, 1>)> f;
        double                                                definite_integral;
    };

    std::vector<TestCase> tests =
    {
        { linspace<double>(0, 2*M_PI, ARR_LEN), [](auto x){ return sin(x); }, 0. },
        { linspace<double>(0., 1., ARR_LEN), [](auto x){ return 3 * x * x; }, 1. },
        { linspace<double>(-1., 1., ARR_LEN), [](auto x){ return 1/(1 + 25*(x*x)); }, 0.54936 },
    };

    for (const auto& test : tests)
    {
        xtensor<double, 1> x = test.x;
        xtensor<double, 1> y = test.f(test.x);
        auto tck = interpolate::splrep(x, y);

        EXPECT_TRUE(allclose(interpolate::splint(test.x[0], test.x[test.x.shape()[0] - 1], tck),
                             test.definite_integral,
                             pow10(-5)));
    }
}

TEST(smoke_test, splrep_spalde)
{
    xtensor<double, 1> x = linspace<double>(0, ARR_LEN/10, ARR_LEN);
    xtensor<double, 1> y = x;

    auto tck = interpolate::splrep(x, y);

    xtensor<double, 1> dx = arange<double>({ ARR_LEN/10 });
    auto dy = interpolate::spalde(dx, tck);

    for (std::size_t i = 0; i < ARR_LEN/10; ++i)
    {
        xtensor<double, 1> expected = { dx[i], 1, 0, 0 };

        EXPECT_TRUE(allclose(expected, flatten(view(dy, range(i, i+1), all()))));
    }
}

}  // namespace xt
