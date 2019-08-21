/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */
#pragma once

#include "MaterialLib/MPL/Property.h"
#include "MaterialLib/MPL/VariableType.h"

namespace MaterialPropertyLib
{

struct IndependentVariable
{
    Variable type;
    VariableType reference_condition; // scalar or vector
    VariableType slope; // scalar or matrix
};

/// The linear property class. This property calculates the linear relationship.
/// The current implementation accepts only the double datatype defined in
/// PropertyDataType.
class LinearProperty final : public Property
{
public:
    /// This constructor accepts single values of double data type defined in
    /// the PropertyDataType definition and sets the protected attribute _value
    /// of the base class Property to that value.
    LinearProperty(PropertyDataType const& property_reference_value,
                   IndependentVariable const& v);

    /// This method computes the value of a property depending linearly on
    /// the value of the given primary variable.
    PropertyDataType value(VariableArray const& variable_array,
                           ParameterLib::SpatialPosition const& /*pos*/,
                           double const /*t*/) const override;

    /// This method will compute the derivative of a property with respect to
    /// the given primary variable.
    PropertyDataType dValue(VariableArray const& variable_array,
                            Variable const primary_variable,
                           ParameterLib::SpatialPosition const& /*pos*/,
                           double const /*t*/) const override;
                           
    /// This method will compute the second derivative of a
    /// property with respect to the given primary variables pv1 and pv2.
    PropertyDataType d2Value(VariableArray const& variable_array,
                             Variable const pv1,
                             Variable const pv2) const override;

private:
    IndependentVariable const _independent_variable;
};
}  // namespace MaterialPropertyLib
