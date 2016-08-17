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
struct ConstantParameter final : public Parameter<T>
{
    /// Construction with single value.
    explicit ConstantParameter(T const& value) : _values({value}) {}

    /// Construction with a tuple.
    /// The given tuple must be non-empty.
    explicit ConstantParameter(std::vector<T> const& values) : _values(values)
    {
        assert(!values.empty());
    }

    bool isTimeDependent() const override { return false; }

    unsigned getNumberOfComponents() const override
    {
        return static_cast<unsigned>(_values.size());
    }

    std::vector<T> const& operator()(
        double const /*t*/, SpatialPosition const& /*pos*/) const override
    {
        return _values;
    }

private:
    std::vector<T> const _values;
};

std::unique_ptr<ParameterBase> createConstantParameter(
    BaseLib::ConfigTree const& config);

}  // ProcessLib

#endif  // PROCESSLIB_CONSTANTPARAMETER_H
