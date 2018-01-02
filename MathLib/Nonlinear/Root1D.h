/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
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
#include <logog/include/logog.hpp>
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

}


//! Use the regula falsi method to find the root of some scalar function of one
//! variable.
template<typename SubType, typename Function>
class RegulaFalsi
{
public:
    //! Initializes finding a root of the \c Function \c f in the interval
    //! [\c a, \c b].
    RegulaFalsi(Function const& f, double a, double b)
        : _f(f), _a(a), _b(b), _fa(f(a)), _fb(f(b))
    {
        static_assert(
                std::is_same<double, decltype(f(0.0))>::value,
                "Using this class for functions that do not return double"
                " involves a lot of casts. Hence it is disabled.");

        if (detail::almost_zero(_fa)) {
            _b = _a;
        } else if (detail::almost_zero(_fb)) {
            _a = _b;
        } else if (detail::same_sign(_fa, _fb)) {
            OGS_FATAL("Regula falsi cannot be done, because the function values"
                " at the interval ends have the same sign.");
        }
    }

    //! Do \c num_steps iteration of regula falsi.
    void step(const unsigned num_steps)
    {
        for (unsigned i=0; i<num_steps; ++i)
        {
            if (_a == _b) return;

            const double s = (_fb - _fa)/(_b - _a);
            const double c = _a - _fa/s;
            const double fc = _f(c);

            if (detail::almost_zero(fc)) {
                _a = _b = c;
                return;
            }
            if (!detail::same_sign(fc, _fb))
            {
                _a = _b; _fa = _fb;
                _b =  c; _fb =  fc;
            } else {
                const double m = SubType::get_m(_fa, _fb, fc);
                _fa *= m;
                _b = c;
                _fb = fc;
            }
        }
    }

    //! Returns the current estimate of the root.
    double getResult() const
    {
        if (_a == _b) return _a;

        const double s = (_fb - _fa)/(_b - _a);
        const double c = _a - _fa/s;

        return c;
    }

    //! Returns the size of the current search interval.
    double getRange() const { return std::fabs(_a - _b); }

private:
    Function const& _f;
    double _a, _b, _fa, _fb;
};


/*! Creates a new instance of \c RegulaFalsi that uses the given modified iteration
 *  method \c SubType.
 *
 * \see https://en.wikipedia.org/wiki/False_position_method#Improvements_in_regula_falsi
 */
template<typename SubType, typename Function>
RegulaFalsi<SubType, Function>
makeRegulaFalsi(Function const& f, double const a, double const b)
{
    return RegulaFalsi<SubType, Function>(f, a, b);
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

}

} // namespace MathLib
