/**
 * \author Norihiro Watanabe
 * \date   2013-08-13
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <array>
#include <numeric>
#include <utility>

namespace MathLib
{

/**
 * Linear function
 * \f[
 *     f(x_1,...,x_k)=a_0+a_1*x_1+...+a_k*x_k
 * \f]
 *
 * \tparam T_TYPE  value type
 * \tparam N_VARS  the number of variables (k)
 */
template <typename T_TYPE, unsigned N_VARS>
class LinearFunction
{
public:
    /**
     * Constructor
     * \param coefficients  an array of coefficients of a linear function.
     * The size of the coefficient array should equal to the number of variables + 1
     */
    explicit LinearFunction(std::array<T_TYPE, N_VARS + 1> coefficients)
        : _coefficients(std::move(coefficients))
    {}

    /**
     * evaluate the function
     * \param x  an array of variables. the size of the array should equal to the number of variables
     */
    T_TYPE operator()(T_TYPE const * const x) const
    {
        return std::inner_product(_coefficients.cbegin()+1, _coefficients.cend(), x, _coefficients.front());
    }

private:
    /// Coefficients of a linear function
    const std::array<T_TYPE, N_VARS+1> _coefficients;
};

} // MathLib
