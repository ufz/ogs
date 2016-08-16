/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESSLIB_CONSTANTPARAMETER_H
#define PROCESSLIB_CONSTANTPARAMETER_H

#include "Parameter.h"

namespace ProcessLib
{
/// Single, constant value parameter.
template <typename T>
struct ConstantParameter final : public Parameter<T> {
    ConstantParameter(T const& value) : _value{{value}} {}
    std::vector<T> const& getTuple(
        double const /*t*/, SpatialPosition const& /*pos*/) const override
    {
        return _value;
    }

private:
    std::vector<T> _value;
};

std::unique_ptr<ParameterBase> createConstantParameter(
    BaseLib::ConfigTree const& config);

}  // ProcessLib

#endif  // PROCESSLIB_CONSTANTPARAMETER_H
