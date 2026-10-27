// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "MaterialLib/MPL/Property.h"
#include "MaterialLib/MPL/VariableType.h"

namespace MathLib
{
class Curve;
}  // namespace MathLib

namespace MaterialPropertyLib
{
using StringOrVariable = std::variant<std::string, Variable>;
/// Definition of a curve property, i.e., a property interpolating pairs of
/// supporting point and value. The pairs and the interpolation are specified
/// by a curve.
/// The current implementation accepts only the double datatype defined in
/// PropertyDataType.
class Curve final : public Property
{
public:
    /// This constructor allows to specify the independent variable of the
    /// curve and the curve itself.
    Curve(std::string name, StringOrVariable const independent_variable,
          MathLib::Curve const& curve);

    /// This method computes the value of a property depending on the value
    /// of the given primary variable via the curve.
    PropertyDataType value(VariableArray const& variable_array,
                           ParameterLib::SpatialPosition const& pos,
                           double const t, double const dt) const override;

    /// This method will compute the derivative of a property with respect to
    /// the given primary variable.
    PropertyDataType dValue(VariableArray const& variable_array,
                            Variable const variable,
                            ParameterLib::SpatialPosition const& /*pos*/,
                            double const /*t*/,
                            double const /*dt*/) const override;

private:
    /// The variable type that the curve property depends on.
    StringOrVariable const independent_variable_;
    /// The curve used by the property.
    MathLib::Curve const& curve_;
};
}  // namespace MaterialPropertyLib
