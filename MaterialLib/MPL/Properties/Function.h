// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <map>
#include <utility>
#include <vector>

#include "MaterialLib/MPL/Property.h"
#include "MaterialLib/MPL/VariableType.h"
#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"

namespace MaterialPropertyLib
{
/// A function property defined by mathematical expression. For the evaluation
/// of the expressions the exprtk library is used. In the expressions all
/// variables defined in MaterialPropertyLib::Variable enum, t for time,
/// x,y,z for the spatial position are supported, and curves from the
/// `<curves>` section can be called using their names. A curve is a single
/// argument function and can be used in an expression like `curveA(sin(t))`.
///
/// The evaluation is thread-safe for OpenMP by using per-thread storage.
class Function final : public Property
{
public:
    Function(
        std::string name,
        std::vector<std::string> const& value_string_expressions,
        std::vector<std::pair<std::string, std::vector<std::string>>> const&
            dvalue_string_expressions,
        std::map<std::string,
                 std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
            curves);

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
};
}  // namespace MaterialPropertyLib
