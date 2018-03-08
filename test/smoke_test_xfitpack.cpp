// xtensor-fitpack: https://github.com/rjsberry/xtensor-fitpack
//
// Copyright (C) 2018, Richard Berry <rjsberry@protonmail.com>
//
// Distributed under the terms of BSD 2-Clause "simplified" license. (See
// accompanying file LICENSE, or copy at
// https://github.com/rjsberry/xtensor-fitpack/blob/master/LICENSE)
//

// Smoke tests for fitpack routines: generally just check they run safely.

#include "gtest/gtest.h"
#include "xtensor/xtensor.hpp"
#include "xtensor/xrandom.hpp"

#include "xtensor-fitpack/xfitpack.hpp"

const std::size_t TEST_ARRAY_LENGTH = 100;

TEST(xfitpack_smoke_test, splrep) {
    xt::xtensor<double, 1> x = xt::arange<double>(TEST_ARRAY_LENGTH);
    xt::xtensor<double, 1> y = xt::random::randn<double>({TEST_ARRAY_LENGTH});

    auto tck = xt::fitpack::splrep(x, y);
}

TEST(xfitpack_smoke_test, splev) {
    xt::xtensor<double, 1> x = xt::arange<double>(TEST_ARRAY_LENGTH);
    xt::xtensor<double, 1> y = xt::random::randn<double>({TEST_ARRAY_LENGTH});

    auto tck = xt::fitpack::splrep(x, y);

    xt::xtensor<double, 1> x_interp = xt::arange<double>(1, TEST_ARRAY_LENGTH-1) - 0.5;
    auto y_interp = xt::fitpack::splev(x_interp, tck);
}
