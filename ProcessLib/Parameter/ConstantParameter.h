/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <utility>

#include "Parameter.h"

namespace ProcessLib
{
/// Single, constant value parameter.
template <typename T>
struct ConstantParameter final : public Parameter<T>
{
    /// Construction with single value.
    explicit ConstantParameter(std::string const& name_, T const& value)
        : Parameter<T>(name_), _values({value})
    {
    }

    /// Construction with a tuple.
    /// The given tuple must be non-empty.
    explicit ConstantParameter(std::string const& name_, std::vector<T> values)
        : Parameter<T>(name_), _values(std::move(values))
    {
        assert(!_values.empty());
    }

    bool isTimeDependent() const override { return false; }

    int getNumberOfComponents() const override
    {
        return static_cast<int>(_values.size());
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
    std::string const& name, BaseLib::ConfigTree const& config);

}  // ProcessLib
