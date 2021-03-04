/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on March 4, 2021, 3:05 PM
 */

#pragma once

#include "MaterialLib/MPL/Property.h"

namespace MaterialPropertyLib
{
class Phase;
class WaterVapourDensity final : public Property
{
public:
    explicit WaterVapourDensity(std::string name) { name_ = std::move(name); }
    void checkScale() const override
    {
        if (!std::holds_alternative<Phase*>(scale_))
        {
            OGS_FATAL(
                "The property 'WaterVapourDensity' is "
                "implemented on the 'Phase' scale only.");
        }
    }

    /// \return The water vapour density.
    PropertyDataType value(VariableArray const& variable_array,
                           ParameterLib::SpatialPosition const& pos,
                           double const t, double const dt) const override;
    /// \return The derivative of  water vapour density with respect to
    /// temperature or phase (water) pressure.
    PropertyDataType dValue(VariableArray const& variable_array,
                            Variable const primary_variable,
                            ParameterLib::SpatialPosition const& pos,
                            double const t, double const dt) const override;
};
}  // namespace MaterialPropertyLib
