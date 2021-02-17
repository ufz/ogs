/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on February 17, 2021, 3:27 PM
 */

#pragma once

#include "BaseLib/Error.h"
#include "MaterialLib/MPL/Property.h"

namespace MaterialPropertyLib
{
class Medium;
class SoilThermalConductivitySomerton final : public Property
{
public:
    SoilThermalConductivitySomerton(std::string name,
                                    double const dry_thermal_conductivity,
                                    double const wet_thermal_conductivity)
        : dry_thermal_conductivity_(dry_thermal_conductivity),
          wet_thermal_conductivity_(wet_thermal_conductivity)
    {
        name_ = std::move(name);
        if (dry_thermal_conductivity_ > wet_thermal_conductivity_)
        {
            OGS_FATAL(
                "In 'SoilThermalConductivitySomerton', "
                "dry_thermal_conductivity of '{:g}' is "
                "larger than wet_thermal_conductivity of '{:g}'.",
                dry_thermal_conductivity_, wet_thermal_conductivity_);
        }
    }

    void checkScale() const override
    {
        if (!std::holds_alternative<Medium*>(scale_))
        {
            OGS_FATAL(
                "The property 'SoilThermalConductivitySomerton' is "
                "implemented on the 'media' scale only.");
        }
    }

    PropertyDataType value(VariableArray const& variable_array,
                           ParameterLib::SpatialPosition const& pos,
                           double const t,
                           double const dt) const override;

    PropertyDataType dValue(VariableArray const& variable_array,
                            Variable const variable,
                            ParameterLib::SpatialPosition const& pos,
                            double const t,
                            double const dt) const override;

private:
    double const dry_thermal_conductivity_;
    double const wet_thermal_conductivity_;
};

}  // namespace MaterialPropertyLib
