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

namespace MaterialPropertyLib
{
struct ExponentData
{
    Variable type;
    VariableType
        reference_condition;  ///< a reference condition value of independent
                              ///< variable, for instance reference temperature.
    VariableType factor;      ///< a dimensionless exponent.
};

/// An exponential property.
///
/// This property calculates the exponential relationship \f$ \alpha(\beta) =
/// \alpha_{\mathrm{offset}} + \alpha_{\mathrm{ref}} \cdot \exp (m (\beta -
/// \beta_{\mathrm{ref}})\f$, where:
///  - \f$\alpha_{\mathrm{ref}}\f$ is a reference value, for instance reference
///    viscosity,
///  - \f$\alpha_{\mathrm{offset}}\f$ is additive offset in units of the
///    property,
///  - \f$m\f$ is a dimensionless exponent,
///  - \f$\beta_{\mathrm{ref}}\f$ is a reference condition value of independent
///  variable, for instance reference temperature.
///
/// The current implementation accepts only the scalar independent variables.
class Exponential final : public Property
{
public:
    /// This constructor accepts single values of double data type defined in
    /// the PropertyDataType definition and sets the protected attribute value_
    /// of the base class Property to that value.
    Exponential(std::string name,
                double const offset,
                PropertyDataType const& property_reference_value,
                ExponentData const& v);
    /// This method computes the value of a property \f$\alpha\f$ depending
    /// exponentialy on the value of the given primary variable \f$\beta\f$.
    PropertyDataType value(VariableArray const& variable_array,
                           ParameterLib::SpatialPosition const& pos,
                           double const t, double const dt) const override;
    /// This method will compute the derivative of a property with respect
    /// to the given primary variable.
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
    ExponentData const exponent_data_;
    double const offset_;  ///< additive offset in units of the property.
};
}  // namespace MaterialPropertyLib
