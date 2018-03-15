// xtensor-interpolate: https://github.com/rjsberry/xtensor-interpolate
//
// Copyright (C) 2018, Richard Berry <rjsberry@protonmail.com>
//
// Distributed under the terms of BSD 2-Clause "simplified" license. (See
// accompanying file LICENSE, or copy at
// https://github.com/rjsberry/xtensor-interpolate/blob/master/LICENSE)
//

#include <string>

#include "gtest/gtest.h"

#define XTENSOR_ENABLE_ASSERT
#include "xtensor/xexception.hpp"
#include "xtensor/xmath.hpp"
#include "xtensor/xrandom.hpp"
#include "xtensor/xtensor.hpp"

#include "xtensor-interpolate/xinterpolate.hpp"

namespace xt
{

TEST(univariate_spline, smoke_test)
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

        auto s = interpolate::UnivariateSpline(x, y).set_smoothing_factor(0.0);

        xtensor<double, 1> x_interp = linspace<double>(-M_PI, M_PI, test_arr_len);
        auto y_interp = s(x_interp);

        EXPECT_TRUE(allclose(y_interp, f(x_interp), 0.005));
    }
}

namespace
{

const std::size_t arr_len = 10;

xtensor<double, 1> x = arange<double>(arr_len);
xtensor<double, 1> y = random::randn<double>({arr_len});

}  // namespace

TEST(univariate_spline, invalid_instantiation)
{
    auto y_trim = view(y, range(0, y.shape()[0] - 1));
    try {
        auto s = interpolate::UnivariateSpline(x, y_trim);
        FAIL() << "Expected `std::runtime_error`";
    } catch (const std::runtime_error& e) {
        EXPECT_NE(std::string(e.what()).find("x_.shape()[0] == y_.shape()[0]"), std::string::npos);
    } catch (...) {
        FAIL() << "Expected `std::runtime_error`";
    }
}

TEST(univariate_spline, set_order)
{
    auto s_k1 = interpolate::UnivariateSpline(x, y).set_order(1);
    auto s_k5 = interpolate::UnivariateSpline(x, y).set_order(5);

    try {
        auto s_k0 = interpolate::UnivariateSpline(x, y).set_order(0);
        FAIL() << "Expected `std::out_of_range`";
    } catch (const std::out_of_range& e) {
        EXPECT_EQ(e.what(), std::string("order must be in {1, 2, 3, 4, 5}"));
    } catch (...) {
        FAIL() << "Expected `std::out_of_range`";
    }

    try {
        auto s_k6 = interpolate::UnivariateSpline(x, y).set_order(6);
        FAIL() << "Expected `std::out_of_range`";
    } catch (const std::out_of_range& e) {
        EXPECT_EQ(e.what(), std::string("order must be in {1, 2, 3, 4, 5}"));
    } catch (...) {
        FAIL() << "Expected `std::out_of_range`";
    }
}

TEST(univariate_spline, set_weights)
{
    xtensor<double, 1> w = random::randn<double>({arr_len});
    auto s_pass = interpolate::UnivariateSpline(x, y).set_weights(w);

    auto w_trim = view(w, range(0, w.shape()[0] - 1));
    try {
       auto s_fail = interpolate::UnivariateSpline(x, y).set_weights(w_trim);
       FAIL() << "Expected `std::runtime_error";
    } catch (const std::runtime_error& e) {
        EXPECT_NE(std::string(e.what()).find("w.shape()[0] == m_"), std::string::npos);
    } catch (...) {
       FAIL() << "Expected `std::runtime_error";
    }
}

TEST(univariate_spline, set_extrapolation_mode_int)
{
    auto s_ext0 = interpolate::UnivariateSpline(x, y).set_extrapolation_mode(0);
    auto s_ext3 = interpolate::UnivariateSpline(x, y).set_extrapolation_mode(3);

    try {
        auto s_ext_below = interpolate::UnivariateSpline(x, y).set_extrapolation_mode(-1);
        FAIL() << "Expected `std::out_of_range`";
    } catch (const std::out_of_range& e) {
        EXPECT_EQ(e.what(), std::string("extrapolation mode must be in {0, 1, 2, 3}"));
    } catch (...) {
        FAIL() << "Expected `std::out_of_range`";
    }

    try {
        auto s_ext_above = interpolate::UnivariateSpline(x, y).set_extrapolation_mode(4);
        FAIL() << "Expected `std::out_of_range`";
    } catch (const std::out_of_range& e) {
        EXPECT_EQ(e.what(), std::string("extrapolation mode must be in {0, 1, 2, 3}"));
    } catch (...) {
        FAIL() << "Expected `std::out_of_range`";
    }
}

TEST(univariate_spline, set_extrapolation_mode_string)
{
    auto s_ext_extrapolate = interpolate::UnivariateSpline(x, y).set_extrapolation_mode("extrapolate");
    auto s_ext_zeros = interpolate::UnivariateSpline(x, y).set_extrapolation_mode("zeros");
    auto s_ext_raise = interpolate::UnivariateSpline(x, y).set_extrapolation_mode("raise");
    auto s_ext_const = interpolate::UnivariateSpline(x, y).set_extrapolation_mode("const");

    try {
        auto s_ext_invalid = interpolate::UnivariateSpline(x, y).set_extrapolation_mode("random");
        FAIL() << "Expected `map::at`";
    } catch (const std::out_of_range& e) {
        EXPECT_EQ(e.what(), std::string("map::at"));
    } catch (...) {
        FAIL() << "Expected `map::at`";
    }
}

TEST(univariate_spline, evaluate_invalid_derivative)
{
    auto order = random::randint<int>({1}, 2, 6)[0];
    auto s = interpolate::UnivariateSpline(x, y).set_order(order);

    xtensor<double, 1> x_eval = view(x, range(1, x.shape()[0] - 1)) + 0.5;
    auto y_nu_below = s(x_eval, order - 1);
    auto y_nu_equal = s(x_eval, order);

    try {
        auto y_nu_above = s(x_eval, order + 1);
        FAIL() << "Expected `std::out_of_range`";
    } catch (const std::out_of_range& e) {
        EXPECT_EQ(e.what(), std::string("derivative must be <= order"));
    } catch (...) {
        FAIL() << "Expected `std::out_of_range`";
    }
}

}  // namespace xt
