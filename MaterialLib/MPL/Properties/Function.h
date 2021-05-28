/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */
#pragma once

#include <exprtk.hpp>
#include <utility>
#include <vector>

#include "MaterialLib/MPL/Property.h"
#include "MaterialLib/MPL/VariableType.h"

namespace MaterialPropertyLib
{
/// A function property defined by mathematical expression. For the evaluation
/// of the expressions the exprtk library is used. In the expressions all
/// variables defined in MaterialPropertyLib::Variable enum can be used.
///
/// \warning The evaluation calls are not to be used in parallel (openMP),
/// because the values' updates are using the same space.
class Function final : public Property
{
public:
    Function(std::string name,
             std::vector<std::string>
                 value_string_expressions,
             std::vector<std::pair<std::string, std::vector<std::string>>>
                 dvalue_string_expressions);

    PropertyDataType value(VariableArray const& variable_array,
                           ParameterLib::SpatialPosition const& pos,
                           double const t,
                           double const dt) const override;

    PropertyDataType dValue(VariableArray const& variable_array,
                            Variable const primary_variable,
                            ParameterLib::SpatialPosition const& pos,
                            double const t,
                            double const dt) const override;

private:
    using Expression = exprtk::expression<double>;

    /// Mapping from variable array index to symbol table values.
    std::vector<std::pair<int, double*>> symbol_values_;
    /// Value expressions.
    /// Multiple expressions are representing vector-valued functions.
    std::vector<Expression> value_expressions_;
    /// Derivative expressions with respect to the variable.
    /// Multiple expressions are representing vector-valued functions.
    std::vector<std::pair<Variable, std::vector<Expression>>>
        dvalue_expressions_;
};
}  // namespace MaterialPropertyLib
