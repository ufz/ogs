/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on Feb 8, 2023, 3:05 PM
 */

#pragma once

#include "MaterialLib/MPL/Property.h"

namespace MaterialPropertyLib
{
class Phase;

/// Water saturation temperature model
/// based on the IAPWS Industrial Formulation 1997
/// <a href="http://www.iapws.org/relguide/IF97-Rev.pdf">IF97-Rev</a>
struct WaterSaturationTemperatureIAPWSIF97Region4 final : public Property
{
    explicit WaterSaturationTemperatureIAPWSIF97Region4(std::string name)
    {
        name_ = std::move(name);
    }
    void checkScale() const override
    {
        if (!std::holds_alternative<Phase*>(scale_))
        {
            OGS_FATAL(
                "The property 'WaterSaturationTemperatureIAPWSIF97Region4' is "
                "implemented on the 'Phase' scale only.");
        }
    }

    /// \return The water saturation temperature.
    PropertyDataType value(VariableArray const& variable_array,
                           ParameterLib::SpatialPosition const& pos,
                           double const t, double const dt) const override;
    /// \return The derivative of water saturation temperature.
    PropertyDataType dValue(VariableArray const& variable_array,
                            Variable const variable,
                            ParameterLib::SpatialPosition const& pos,
                            double const t, double const dt) const override;
};
}  // namespace MaterialPropertyLib
