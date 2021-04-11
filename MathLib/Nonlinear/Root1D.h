/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <cassert>
#include <cmath>
#include <limits>
#include <type_traits>
#include <utility>

#include "BaseLib/Error.h"

namespace MathLib
{
namespace Nonlinear
{

namespace detail
{

//! Tells if \c a and \c b have the same sign.
inline bool same_sign(double a, double b)
{
    // note: signbit() casts integers to double and distinguishes between +0.0 and
    // -0.0. However, the latter is not a problem, because if f(x0) == 0, we've
    // already found the root x0 that we are searching for.
    return std::signbit(a) == std::signbit(b);
}

inline bool almost_zero(double a)
{
    return std::abs(a) <= std::numeric_limits<double>::epsilon();
}

}  // namespace detail

//! Use the regula falsi method to find the root of some scalar function of one
//! variable.
template<typename SubType, typename Function>
class RegulaFalsi
{
public:
    //! Initializes finding a root of the \c Function \c f in the interval
    //! [\c a, \c b].
    RegulaFalsi(Function&& f, double a, double b)
        : f_(f), a_(a), b_(b), fa_(f(a)), fb_(f(b))
    {
        static_assert(std::is_same_v<double, decltype(f(0.0))>,
                      "Using this class for functions that do not return double"
                      " involves a lot of casts. Hence it is disabled.");

        if (detail::almost_zero(fa_))
        {
            b_ = a_;
        }
        else if (detail::almost_zero(fb_))
        {
            a_ = b_;
        }
        else if (detail::same_sign(fa_, fb_))
        {
            OGS_FATAL("Regula falsi cannot be done, because the function values"
                " at the interval ends have the same sign.");
        }
    }

    //! Do \c num_steps iteration of regula falsi.
    void step(const unsigned num_steps)
    {
        for (unsigned i=0; i<num_steps; ++i)
        {
            if (a_ == b_)
            {
                return;
            }

            const double s = (fb_ - fa_) / (b_ - a_);
            const double c = a_ - fa_ / s;
            const double fc = f_(c);

            if (detail::almost_zero(fc)) {
                a_ = b_ = c;
                return;
            }
            if (!detail::same_sign(fc, fb_))
            {
                a_ = b_;
                fa_ = fb_;
                b_ = c;
                fb_ = fc;
            } else {
                const double m = SubType::get_m(fa_, fb_, fc);
                fa_ *= m;
                b_ = c;
                fb_ = fc;
            }
        }
    }

    //! Returns the current estimate of the root.
    double getResult() const
    {
        if (a_ == b_)
        {
            return a_;
        }

        const double s = (fb_ - fa_) / (b_ - a_);
        const double c = a_ - fa_ / s;

        return c;
    }

    //! Returns the size of the current search interval.
    double getRange() const { return std::abs(a_ - b_); }

private:
    Function f_;
    double a_, b_, fa_, fb_;
};


/*! Creates a new instance of \c RegulaFalsi that uses the given modified iteration
 *  method \c SubType.
 *
 * \see https://en.wikipedia.org/wiki/False_position_method#Improvements_in_regula_falsi
 */
template <typename SubType, typename Function>
RegulaFalsi<SubType, Function> makeRegulaFalsi(Function&& f,
                                               double const a,
                                               double const b)
{
    return RegulaFalsi<SubType, Function>(std::forward<Function>(f), a, b);
}


//! Used by \c RegulaFalsi in the original regula falsi algorithm.
struct Unmodified
{
    static double get_m(const double /*fa*/, const double /*fb*/, const double /*fc*/)
    { return 1.0; }
};

//! Used by \c RegulaFalsi in a modified version of the regula falsi algorithm.
struct Illinois
{
    static double get_m(const double /*fa*/, const double /*fb*/, const double /*fc*/)
    { return 0.5; }
};

//! Used by \c RegulaFalsi in a modified version of the regula falsi algorithm.
struct Pegasus
{
    static double get_m(const double /*fa*/, const double fb, const double fc)
    { return fb / (fb+fc); }
};

//! Used by \c RegulaFalsi in a modified version of the regula falsi algorithm.
struct AndersonBjorck
{
    static double get_m(const double /*fa*/, const double fb, const double fc)
    {
        const double f = 1.0 - fc / fb;
        return (f >= 0.0) ? f : 0.5;
    }
};

}  // namespace Nonlinear

} // namespace MathLib
