/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on March 19, 2021, 11:49 AM
 */

#pragma once

#include "MaterialLib/MPL/Property.h"
#include "MaterialLib/PhysicalConstant.h"

namespace MaterialPropertyLib
{
class Phase;

class GeneralWaterVapourLatentHeat final : public Property
{
public:
    explicit GeneralWaterVapourLatentHeat(std::string name)
    {
        name_ = std::move(name);
    }

    void checkScale() const override
    {
        if (!std::holds_alternative<Phase*>(scale_))
        {
            OGS_FATAL(
                "The property 'GeneralWaterVapourLatentHeat' is "
                "implemented on the 'phase' scale only.");
        }
    }

    PropertyDataType value(VariableArray const& variable_array,
                           ParameterLib::SpatialPosition const& pos,
                           double const t,
                           double const dt) const override;

    PropertyDataType dValue(VariableArray const& variable_array,
                            Variable const primary_variable,
                            ParameterLib::SpatialPosition const& pos,
                            double const t, double const dt) const override;

private:
    /// Critical temperature.
    static double constexpr T_c_ =
        373.92 + MaterialLib::PhysicalConstant::CelsiusZeroInKelvin;
    ;
};

}  // namespace MaterialPropertyLib
