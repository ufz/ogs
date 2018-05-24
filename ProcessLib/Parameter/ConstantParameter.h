/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
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

    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> getNodalValuesOnElement(
        MeshLib::Element const& element, double const /*t*/) const override
    {
        auto const n_nodes = element.getNumberOfNodes();
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> result(
            n_nodes, getNumberOfComponents());

        // Column vector of values, copied for each node.
        auto const row_values =
            Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, 1> const>(
                _values.data(), _values.size());
        for (unsigned i = 0; i < n_nodes; ++i)
        {
            result.row(i) = row_values;
        }
        return result;
    }

private:
    std::vector<T> const _values;
};

std::unique_ptr<ParameterBase> createConstantParameter(
    std::string const& name, BaseLib::ConfigTree const& config);

}  // ProcessLib
