#ifndef INCLUDE_XTENSOR_INTERPOLATE_XSPLINE_HPP_
#define INCLUDE_XTENSOR_INTERPOLATE_XSPLINE_HPP_

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
#include "xtensor-interpolate/xinterpolate.hpp"

namespace xt
{
namespace interpolate
{

namespace detail
{

auto fpcurf0(const xtensor<double, 1>& x, const xtensor<double, 1>& y,
             const xtensor<double, 1>& w, int m, double xb, double xe, int k,
             double s, int nest, int& n, xtensor<double, 1>& t,
             xtensor<double, 1>& c, double& fp, xtensor<double, 1>& fpint,
             xtensor<int, 1>& nrdata, int iopt = 0, double tol = 0.001,
             int maxit = 20)
{
    auto k1 = k + 1;
    auto k2 = k + 2;
    xtensor<double, 1> wrk = zeros<double>({ nest*3*k2 + m*k1 });
    xtensor<double, 1> wrkn = wrk + nest;
    xtensor<double, 1> wrknk2 = wrk + nest*k2;
    xtensor<double, 1> wrkn2k2 = wrk + nest*2*k2;
    xtensor<double, 1> wrkn3k2 = wrk + nest*3*k2;

    auto ier = 0;

    fp_fpcurf(&iopt, &x[0], &y[0], &w[0], &m, &xb, &xe, &k, &s, &nest, &tol,
              &maxit, &k1, &k2, &n, &t[0], &c[0], &fp, &fpint[0], &wrk[0],
              &wrkn[0], &wrknk2[0], &wrkn2k2[0], &wrkn3k2[0], &nrdata[0], &ier);

    return ier;
}

}  // namespace detail

class Spline
{
  public:
    virtual inline xtensor<double, 1> operator() (const xtensor<double, 1>& x,
                                                  int nu = 0) = 0;
};

/// One-dimensional smoothing spline fit to a given set of data points.
///
/// **Example usage:**
/// @code
/// #include <cmath>
/// #include "xtensor/xtensor.hpp"
/// #include "xtensor/xmath.hpp"
/// #include "xtensor-interpolate/xinterpolate.hpp"
///
/// int main(void) {
///     xt::xtensor<double, 1> x = xt::linspace(0, 2*M_PI, 100);
///     xt::xtensor<double, 1> y = xt::sin(x);
///
///     auto s = xt::interpolate::UnivariateSpline(x, y)
///                                   .set_order(3)
///                                   .set_smoothing_factor(0.0);
///
///     xt::xtensor<double, 1> xs = xt::linspace(0, 2*M_PI, 1000);
///     auto ys = s(xs);
///
///     return 0;
/// }
/// @endcode
///
class UnivariateSpline : public virtual Spline
{
  protected:

    // Required interpolation parameters.
    const xtensor<double, 1> x_;
    const xtensor<double, 1> y_;
    xtensor<double, 1>       w_;
    int                      m_;
    double                   xb_;
    double                   xe_;
    int                      k_;
    double                   s_;
    int                      ext_;

    // Variable data (subject to modification by FITPACK).
    bool               req_eval_;
    int                nest_;
    int                n_;
    xtensor<double, 1> t_;
    xtensor<double, 1> c_;
    double             fp_;
    xtensor<double, 1> fpint_;
    xtensor<int, 1>    nrdata_;

  public:

    /// Construct a `UnivariateSpline` object.
    ///
    /// @param [in] x
    ///     1-D array of input data. Must be increasing.
    /// @param [in] y
    ///     1-D array of input data corresponding to `x`. Must be the same length
    ///     as `x`.
    ///
    /// @throws std::runtime_error
    ///     Thrown if `x.shape()[0] != y.shape()[0]`.
    ///
    inline UnivariateSpline(const xtensor<double, 1>& x, const xtensor<double, 1>& y)
        : x_(x)
        , y_(y)
        , w_(ones<double>({ x.shape()[0] }))
        , m_(static_cast<int>(x.shape()[0]))
        , xb_(x[0])
        , xe_(x[x.shape()[0] - 1])
        , k_(3)
        , s_(static_cast<double>(x.shape()[0]))
        , ext_(0)
        , req_eval_(true)
    {
        XTENSOR_ASSERT(x_.shape()[0] == y_.shape()[0]);
    }

    /// @name Attribute Setters
    ///
    /// Optional setters to be chained during object instantiation.
    ///
    /// Can also be used to alter spline object parameters after
    /// instantiation, however, this will require re-calculation of spline
    /// properties with *FITPACK*.
    ///
    /// @{

    /// Set weights for spline fitting.
    ///
    /// @param [in] w
    ///     1-D array of positive weights. Must have the same length as `x` and
    ///     `y`. Defaults to an array of ones.
    ///
    /// @throws std::runtime_error
    ///     Thrown if `w.shape()[0] != m` where `m` is the length of both `x`
    ///     and `y`.
    ///
    inline auto set_weights(const xtensor<double, 1>& w)
    {
        XTENSOR_ASSERT(w.shape()[0] == m_);
        w_ = w;
        return reset(*this);
    }

    /// Set the boundary of the approximation interval.
    ///
    /// @param [in] begin
    ///     Defaults to `x[0]`.
    /// @param [in] end
    ///     Defaults to `x[x.shape()[0] - 1]`
    ///
    inline auto set_bounds(double begin, double end)
    {
        xb_ = begin;
        xe_ = end;
        return reset(*this);
    }
    
    /// Set the degree of the smoothing spline.
    ///
    /// @param [in] k
    ///     Degree of the smoothing spline. Defaults to `3`.
    ///
    /// @throws std::out_of_range
    ///     Thrown if `k` does not satisfy `1 <= k <= 5`.
    ///
    inline auto set_order(int k)
    {
        if (!(1 <= k && k <= 5))
        {
            throw std::out_of_range("order must be in {1, 2, 3, 4, 5}");
        }
        k_ = k;
        return reset(*this);
    }

    /// Set the positive smoothing factor used to choose the number of knots.
    ///
    /// @param [in] s
    ///     Defaults to `w.shape()[0]`.
    ///     If `0`, spline will interpolate through all data points.
    ///
    /// @throws std::out_of_range
    ///     Thrown if `s < 0`.
    ///
    inline auto set_smoothing_factor(double s)
    {
        s_ = s;
        return reset(*this);
    }

    /// Set the extrapolation mode via integer.
    ///
    /// @param [in] ext
    ///     Defaults to `0`. Options:
    ///       - `0`: Return the extrapolated values.
    ///       - `1`: Return 0.
    ///       - `2`: Throw an exception.
    ///       - `3`: Return the boundary value.
    ///
    /// @throws std::out_of_range
    ///     Thrown if `ext` not recognized.
    ///
    inline auto set_extrapolation_mode(int ext)
    {
        if (!(0 <= ext && ext <= 3))
        {
            throw std::out_of_range("extrapolation mode must be in {0, 1, 2, 3}");
        }
        ext_ = ext;
        return reset(*this);
    }

    /// Set the extrapolation mode via string.
    ///
    /// @param [in] ext
    ///     Defaults to `"extrapolate"`. Options:
    ///       - `"extrapolate"`: Return the extrapolated values.
    ///       - `"zeros"`: Return 0.
    ///       - `"raise"`: Throw an exception.
    ///       - `"const"`: Return the boundary value.
    ///
    /// @throws std::out_of_range
    ///     Thrown if `ext` not recognized.
    ///
    inline auto set_extrapolation_mode(const std::string& ext)
    {
        std::map<std::string, int> ext_map =
        {
            {"extrapolate", 0}, {"zeros", 1}, {"raise", 2}, {"const", 3}
        };
        ext_ = ext_map.at(ext);
        return reset(*this);
    }

    /// @}  // Attribute Setters

    /// Evaluate the spline (or it's `nu`-th derivate) at positions given by x.
    ///
    /// @param [in] x
    ///     An array of points at which to return the value of the smoothed
    ///     spline or its derivatives.
    /// @param [in] nu
    ///     The order of derivative of the spline to compute (must be less than
    ///      or equal to the spline order).
    ///
    /// @throw std::out_of_range
    ///     Thrown if the requested derivative is invalid.
    ///
    inline xtensor<double, 1> operator() (const xtensor<double, 1>& x, int nu = 0)
    {
        if (!(nu <= k_))
        {
            throw std::out_of_range("derivative must be <= order");
        }
        if (!x.shape()[0])
        {
            return {};
        }
        return splev(x, get_tck(), nu, ext_);
    }

    /// Return the definite integral of the spline between two points.
    ///
    /// @param [in] a
    ///     Lower limit of integration.
    /// @param [in] b
    ///     Upper limit of integration.
    ///
    /// @note This routine assumes that the spline is **0** outside of its data limits.
    ///
    inline double integral(double a, double b)
    {
        return splint(a, b, get_tck());
    }

    /// Return all derivatives of the spline at point `x`.
    ///
    /// @param [in] x
    ///     The point at which to evaluate derivates.
    ///
    inline auto derivatives(double x)
    {
        return spalde(x_, get_tck());
    }

  private:

    /// Evaluate `tck` based on spline parameters.
    ///
    void evaluate(void)
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

        req_eval_ = false;
    }

    /// Get the `tck` tuple used for spline interpolation.
    ///
    inline std::tuple<xtensor<double, 1>, xtensor<double, 1>, int> get_tck(void)
    {
        if (req_eval_)
        {
            evaluate();
        }
        return std::make_tuple(t_, c_, k_);
    }

    /// Get a deinitialized version of the class.
    ///
    inline UnivariateSpline& reset(const UnivariateSpline& self)
    {
        req_eval_ = true;
        return *this;
    }
};

}  // namespace interpolate

}  // namespace xt

#endif  // INCLUDE_XTENSOR_INTERPOLATE_XSPLINE_HPP_
