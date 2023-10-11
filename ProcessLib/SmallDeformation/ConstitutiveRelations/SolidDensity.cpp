/**
 * \file
 * \copyright
 * Copyright (c) 2012-2023, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "SolidDensity.h"

namespace ProcessLib::SmallDeformation
{
void SolidDensityModel::eval(SpaceTimeData const& x_t,
                             MediaData const& media_data,
                             TemperatureData<DisplacementDim> const& T_data,
                             SolidDensity& out) const
{
    namespace MPL = MaterialPropertyLib;
    MPL::VariableArray variables;
    variables.temperature = T_data.T;

    *out = media_data.solid.property(MPL::PropertyType::density)
               .template value<double>(variables, x_t.x, x_t.t, x_t.dt);
}
}  // namespace ProcessLib::SmallDeformation
