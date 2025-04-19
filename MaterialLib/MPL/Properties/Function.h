/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */
#pragma once

#include <utility>
#include <vector>

#include "MaterialLib/MPL/Property.h"
#include "MaterialLib/MPL/VariableType.h"

namespace MaterialPropertyLib
{
/// A function property defined by mathematical expression. For the evaluation
/// of the expressions the exprtk library is used. In the expressions all
/// variables defined in MaterialPropertyLib::Variable enum, t for time and
/// x,y,z for the spatial position can be used.
///
/// \warning The evaluation calls are not to be used in parallel (openMP),
/// because the values' updates are using the same space.
class Function final : public Property
{
public:
    Function(
        std::string name,
        std::vector<std::string> const& value_string_expressions,
        std::vector<std::pair<std::string, std::vector<std::string>>> const&
            dvalue_string_expressions);

    PropertyDataType value(VariableArray const& variable_array,
                           ParameterLib::SpatialPosition const& pos,
                           double const t,
                           double const dt) const override;

    PropertyDataType dValue(VariableArray const& variable_array,
                            Variable const variable,
                            ParameterLib::SpatialPosition const& pos,
                            double const t,
                            double const dt) const override;

    ~Function();

private:
    template <int D>
    class Implementation;

    std::unique_ptr<Implementation<2>> impl2_;
    std::unique_ptr<Implementation<3>> impl3_;

    std::variant<Function::Implementation<2>*, Function::Implementation<3>*>
    getImplementationForDimensionOfVariableArray(
        VariableArray const& variable_array) const;

    /// Variables used in the exprtk expressions.
    std::vector<Variable> variables_;

    mutable std::mutex mutex_;
};
}  // namespace MaterialPropertyLib
