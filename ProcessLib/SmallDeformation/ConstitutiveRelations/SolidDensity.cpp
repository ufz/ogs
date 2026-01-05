// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "SolidDensity.h"

namespace ProcessLib::SmallDeformation
{
void SolidDensityModel::eval(SpaceTimeData const& x_t,
                             MediaData const& media_data,
                             Temperature const& temperature,
                             SolidDensity& out) const
{
    namespace MPL = MaterialPropertyLib;
    MPL::VariableArray variables;
    variables.temperature = *temperature;

    *out = media_data.solid.property(MPL::PropertyType::density)
               .template value<double>(variables, x_t.x, x_t.t, x_t.dt);
}
}  // namespace ProcessLib::SmallDeformation
