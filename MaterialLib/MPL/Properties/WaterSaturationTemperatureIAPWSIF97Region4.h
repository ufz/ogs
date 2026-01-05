// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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
