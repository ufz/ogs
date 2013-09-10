/**
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MATHLIB_WEIGHTEDSUM_H_
#define MATHLIB_WEIGHTEDSUM_H_

namespace MathLib
{

namespace detail {

template <unsigned I, typename Method>
struct SUM
{
    template <typename F>
    static
    double
    add(F const& f)
    {
        return f(Method::X[I - 1]) * Method::W[I - 1] + SUM<I - 1, Method>::add(f);
    }
};

template <typename Method>
struct SUM<0, Method>
{
    template <typename F>
    static
    double
    add(F const&)
    {
        return 0;
    }
};

}   // namespace detail

/// Computes weighted sum using given integration method.
/// The weighted sum over all positions x_i in the integration method is computed
/// as follows:
/// \\sum_{i = 0..Method::Order} (f(x_i) * w_i).
///
/// \tparam Method  Integration method.
/// \tparam Func    Function type.
template <typename Method>
struct WeightedSum
{
    template <typename Func>
    static
    double
    add(Func const& f)
    {
        return detail::SUM<Method::Order, Method>::add(f);
    }
};

}   // namespace MathLib

#endif  // MATHLIB_WEIGHTEDSUM_H_
