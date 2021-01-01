/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */
#pragma once

#include "MaterialLib/MPL/Property.h"
#include "MaterialLib/MPL/VariableType.h"

#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"

namespace MaterialPropertyLib
{
/// Definition of a curve property, i.e., a piecewise linear property given by
/// pairs of supporting point and value. The pairs are specified by a curve.
/// The current implementation accepts only the double datatype defined in
/// PropertyDataType.
class Curve final : public Property
{
public:
    /// This constructor allows to specify the independent variable of the
    /// curve and the curve itself, i.e., a piecewise linear function.
    Curve(std::string name, Variable const independent_variable,
          MathLib::PiecewiseLinearInterpolation const& curve);

    /// This method computes the value of a property depending linearly on
    /// the value of the given primary variable.
    PropertyDataType value(VariableArray const& variable_array,
                           ParameterLib::SpatialPosition const& pos,
                           double const t, double const dt) const override;

    /// This method will compute the derivative of a property with respect to
    /// the given primary variable.
    PropertyDataType dValue(VariableArray const& variable_array,
                            Variable const primary_variable,
                            ParameterLib::SpatialPosition const& /*pos*/,
                            double const /*t*/,
                            double const /*dt*/) const override;

private:
    /// The variable type that the curve property depends on.
    Variable const independent_variable_;
    /// The curve used by the property.
    MathLib::PiecewiseLinearInterpolation const& curve_;
};
}  // namespace MaterialPropertyLib
