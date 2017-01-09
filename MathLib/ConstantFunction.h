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

namespace MathLib
{

/**
 * A constant function.
 * \f[
 *     f(x_1,...,x_k)=a_0
 * \f]
 */
template <typename T_TYPE>
class ConstantFunction
{
public:
    explicit ConstantFunction(T_TYPE const& value)
    : _value(value)
    {}

    T_TYPE operator()() const
    {
        return _value;
    }
private:
    T_TYPE const _value;
};

}
