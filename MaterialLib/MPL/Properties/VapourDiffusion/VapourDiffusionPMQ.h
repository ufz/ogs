/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on March 7, 2021, 9:14 AM
 */

#pragma once

#include "MaterialLib/MPL/Property.h"

namespace MaterialPropertyLib
{
class Phase;

class VapourDiffusionPMQ final : public Property
{
public:
    explicit VapourDiffusionPMQ(std::string name) { name_ = std::move(name); }

    void checkScale() const override
    {
        if (!std::holds_alternative<Phase*>(scale_))
        {
            OGS_FATAL(
                "The property 'VapourDiffusionPMQ' is "
                "implemented on the 'phase' scale only.");
        }
    }

    PropertyDataType value(VariableArray const& variable_array,
                           ParameterLib::SpatialPosition const& pos,
                           double const t,
                           double const dt) const override;

    PropertyDataType dValue(VariableArray const& variable_array,
                            Variable const variable,
                            ParameterLib::SpatialPosition const& pos,
                            double const t, double const dt) const override;

private:
    static double constexpr tortuosity_ = 0.66;
};

}  // namespace MaterialPropertyLib
