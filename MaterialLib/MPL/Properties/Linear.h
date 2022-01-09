/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */
#pragma once

#include "MaterialLib/MPL/Property.h"
#include "MaterialLib/MPL/VariableType.h"

namespace MaterialPropertyLib
{
struct IndependentVariable
{
    Variable type;
    VariableType reference_condition;  // scalar or vector
    VariableType slope;                // scalar or matrix
};

/// The linear property class. This property calculates the linear relationship.
/// The current implementation accepts only the double datatype defined in
/// PropertyDataType.
class Linear final : public Property
{
public:
    /// This constructor accepts single values of double data type defined in
    /// the PropertyDataType definition and sets the protected attribute value_
    /// of the base class Property to that value.
    Linear(std::string name,
           PropertyDataType const& property_reference_value,
           std::vector<IndependentVariable> const& vs);

    /// This method computes the value of a property depending linearly on
    /// the value of the given primary variable.
    PropertyDataType value(VariableArray const& variable_array,
                           ParameterLib::SpatialPosition const& /*pos*/,
                           double const /*t*/,
                           double const /*dt*/) const override;

    /// This method will compute the derivative of a property with respect to
    /// the given primary variable.
    PropertyDataType dValue(VariableArray const& variable_array,
                            Variable const primary_variable,
                            ParameterLib::SpatialPosition const& /*pos*/,
                            double const /*t*/,
                            double const /*dt*/) const override;

    /// This method will compute the second derivative of a
    /// property with respect to the given primary variables pv1 and pv2.
    PropertyDataType d2Value(VariableArray const& variable_array,
                             Variable const pv1, Variable const pv2,
                             ParameterLib::SpatialPosition const& /*pos*/,
                             double const /*t*/,
                             double const /*dt*/) const override;

private:
    std::vector<IndependentVariable> const independent_variables_;
};
}  // namespace MaterialPropertyLib
