/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */
#pragma once

#include "MaterialLib/MPL/Property.h"

namespace MaterialPropertyLib
{
/// Saturation dependent model for effective heat conduction
/// \details This model is for the porous media with isotropic thermal conductivity. This property must be a medium property, it
/// computes the effective thermal conductivity based on a wet
/// and a dry value
/// \f$ K_{\mathrm{eff}} = S K_{\mathrm{wet}} + (1-S) K_{\mathrm{dry}} \f$
class SaturationDependentThermalConductivity final : public Property
{
public:
    SaturationDependentThermalConductivity(std::string name,
                                      double const K_dry,
                                      double const K_wet);

    void checkScale() const override;

    PropertyDataType value(VariableArray const& variable_array,
                           ParameterLib::SpatialPosition const& /*pos*/,
                           double const /*t*/,
                           double const /*dt*/) const override;
    PropertyDataType dValue(VariableArray const& variable_array,
                            Variable const variable,
                            ParameterLib::SpatialPosition const& /*pos*/,
                            double const /*t*/,
                            double const /*dt*/) const override;

private:
    double const K_dry_; //< Effective thermal conductivity of the dry material.
    double const K_wet_; //< Effective thermal conductivity of the wet material.
};
}  // namespace MaterialPropertyLib
