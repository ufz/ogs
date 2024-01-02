/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
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

/** Water temperature in IAPWS-IF97 region1 model
 *  based on the IAPWS Industrial Formulation 1997
 *  <a href="http://www.iapws.org/relguide/IF97-Rev.pdf">IF97-Rev</a>
 */
class WaterTemperatureIAPWSIF97Region1 final : public Property
{
public:
    explicit WaterTemperatureIAPWSIF97Region1(std::string name)
    {
        name_ = std::move(name);
    }
    void checkScale() const override
    {
        if (!std::holds_alternative<Phase*>(scale_))
        {
            OGS_FATAL(
                "The property 'WaterTemperatureIAPWSIF97Region1' is "
                "implemented on the 'Phase' scale only.");
        }
    }

    /// \return The water temperature in IAPWS-IF97 region1.
    PropertyDataType value(VariableArray const& variable_array,
                           ParameterLib::SpatialPosition const& pos,
                           double const t, double const dt) const override;
    /// \return The derivative of water temperature in IAPWS-IF97 region1.
    PropertyDataType dValue(VariableArray const& variable_array,
                            Variable const variable,
                            ParameterLib::SpatialPosition const& pos,
                            double const t, double const dt) const override;
};
}  // namespace MaterialPropertyLib
