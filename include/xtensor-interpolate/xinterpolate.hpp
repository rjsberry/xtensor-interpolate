#ifndef INCLUDE_XTENSOR_INTERPOLATE_HPP_
#define INCLUDE_XTENSOR_INTERPOLATE_HPP_

// xtensor-interpolate: https://github.com/rjsberry/xtensor-interpolate
//
// Copyright (C) 2018, Richard Berry <rjsberry@protonmail.com>
//
// Distributed under the terms of BSD 2-Clause "simplified" license. (See
// accompanying file LICENSE, or copy at
// https://github.com/rjsberry/xtensor-interpolate/blob/master/LICENSE)
//

#include <tuple>

#include "xtensor/xtensor.hpp"
#include "xtensor/xview.hpp"

#include "xfitpack/xfcdec.h"

namespace xt
{
namespace interpolate
{

namespace detail
{

auto fpcurf0(xtensor<double, 1>& x, xtensor<double, 1>& y, xtensor<double, 1>& w,
             int& m, double& xb, double& xe, int& k, double& s, int& nest,
             int& n, xtensor<double, 1>& t, xtensor<double, 1>& c, double& fp,
             xtensor<double, 1>& fpint, xtensor<int, 1>& nrdata,
             int iopt = 0, double tol = 0.001, int maxit = 20)
{
    // Internal data allocations for use in Fortran.
    auto k1 = k + 1;
    auto k2 = k + 2;
    xtensor<double, 1> wrk = zeros<double>({ nest*3*k2 + m*k1 });
    xtensor<double, 1> wrkn = wrk + nest;
    xtensor<double, 1> wrknk2 = wrk + nest*k2;
    xtensor<double, 1> wrkn2k2 = wrk + nest*2*k2;
    xtensor<double, 1> wrkn3k2 = wrk + nest*3*k2;

    auto ier = 0;

    _fc_fpcurf(&iopt, &x[0], &y[0], &w[0], &m, &xb, &xe, &k, &s, &nest,
               &tol, &maxit, &k1, &k2, &n, &t[0], &c[0], &fp, &fpint[0],
               &wrk[0], &wrkn[0], &wrknk2[0], &wrkn2k2[0], &wrkn3k2[0],
               &nrdata[0], &ier);

    return ier;
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
template <class E, class... Args>
auto splev(xexpression<E>& x, std::tuple<Args...>& tck, int der = 0, int ext = 0)
{
    auto _x = x.derived_cast();
    auto m = static_cast<int>(x.derived_cast().shape()[0]);

    auto t = std::get<0>(tck);
    auto n = static_cast<int>(t.size());
    auto c = std::get<1>(tck);
    auto k = std::get<2>(tck);

    if (!(0 <= der && der <= k))
    {
        throw std::out_of_range("der");
    }
    if (!(0 <= ext && ext <= 3))
    {
        throw std::out_of_range("ext");
    }

    xtensor<double, 1> y = zeros<double>({ static_cast<std::size_t>(m) });
    auto ier = 0;

    _fc_splev(&t[0], &n, &c[0], &k, &_x[0], &y[0], &m, &ext, &ier);

    return y;
}

}  // namespace detail

class Spline
{
  public:
    virtual xtensor<double, 1> operator() (xtensor<double, 1>& x, int nu = 0) = 0;
};

// One-dimensional smoothing spline fit to a given set of data points.
//
class UnivariateSpline : public Spline
{
  public:

    // Construct a `UnivariateSpline` object.
    //
    // @param [in] x, y
    //     The data points defining a curve y = f(x).
    //
    UnivariateSpline(const xtensor<double, 1>& x, const xtensor<double, 1>& y);

    // Optional setters to be chained during object instantiation.
    //
    // Can also be used to alter spline object parameters after initial
    // instantiation, however, this will require re-calculation of spline
    // properties using the FC interface.
    //
    inline auto set_weights(const xtensor<double, 1>& w);
    inline auto set_bbox(const xtensor<double, 1>& bbox);
    inline auto set_order(int k);
    inline auto set_smoothing_factor(double s);

    // Evaluate the spline (or it's nu-th derivate) at positions given by x.
    //
    xtensor<double, 1> operator() (xtensor<double, 1>& x, int nu = 0);

  private:

    // Re-evaluate `tck` based on spline parameters.
    //
    void initialize(void);

    enum State
    {
        UNINITIALIZED,
        INITIALIZED
    };

    State state_;

    // Required interpolation parameters.
    xtensor<double, 1> x_;
    xtensor<double, 1> y_;
    xtensor<double, 1> w_;
    int                m_;
    double             xb_;
    double             xe_;
    int                k_;
    double             s_;

    // Variable data (subject to modification by FITPACK).
    int                nest_;
    int                n_;
    xtensor<double, 1> t_;
    xtensor<double, 1> c_;
    double             fp_;
    xtensor<double, 1> fpint_;
    xtensor<int, 1> nrdata_;
};

UnivariateSpline::UnivariateSpline(const xtensor<double, 1>& x,
                                   const xtensor<double, 1>& y)
    : x_(x)
    , y_(y)
    , w_(ones<double>({ x.shape()[0] }))
    , m_(static_cast<int>(x.shape()[0]))
    , xb_(x[0])
    , xe_(x[x.shape()[0] - 1])
    , k_(3)
    , s_(static_cast<double>(x.shape()[0]))
    , state_(UNINITIALIZED)
{
    XTENSOR_ASSERT(x_.shape()[0] == y_.shape()[0]);
}

inline auto UnivariateSpline::set_weights(const xtensor<double, 1>& w)
{
    XTENSOR_ASSERT(w.shape()[0] == m_);
    w_ = w;
    state_ = UNINITIALIZED;
    return *this;
}

inline auto UnivariateSpline::set_bbox(const xtensor<double, 1>& bbox)
{
    xb_ = bbox[0];
    xe_ = bbox[1];
    state_ = UNINITIALIZED;
    return *this;
}

inline auto UnivariateSpline::set_order(int k)
{
    k_ = k;
    state_ = UNINITIALIZED;
    return *this;
}

inline auto UnivariateSpline::set_smoothing_factor(double s)
{
    s_ = s;
    state_ = UNINITIALIZED;
    return *this;
}

xtensor<double, 1> UnivariateSpline::operator() (xtensor<double, 1>& x, int nu)
{
    if (!x.shape()[0])
    {
        return {};
    }
    if (state_ == UNINITIALIZED)
    {
        initialize();
    }
    auto tck = std::make_tuple(t_, c_, k_);
    return detail::splev(x, tck, nu);
}

void UnivariateSpline::initialize(void)
{
    nest_ = (s_) ? std::max(m_/2, 2*(k_+1)) : m_ + k_ + 1;
    XTENSOR_ASSERT(nest_ >= 2*(k_+1));

    n_ = 0;
    t_ = zeros<double>({ nest_ });
    c_ = zeros<double>({ nest_ });
    fp_ = 0;
    fpint_ = zeros<double>({ nest_ });
    nrdata_ = zeros<int>({ nest_ });

    auto ier = detail::fpcurf0(x_, y_, w_, m_, xb_, xe_, k_, s_, nest_,
                               n_, t_, c_, fp_, fpint_, nrdata_);

    t_ = view(t_, range(0, n_));
    c_ = view(c_, range(0, n_));

    state_ = INITIALIZED;
}

}  // namespace interpolate

}  // namespace xt

#endif  // INCLUDE_XTENSOR_INTERPOLATE_HPP_
