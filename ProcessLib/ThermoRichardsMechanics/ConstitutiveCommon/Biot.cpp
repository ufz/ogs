/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "Biot.h"

namespace ProcessLib::ThermoRichardsMechanics
{
void BiotModel::eval(SpaceTimeData const& x_t, MediaData const& media_data,
                     BiotData& out) const
{
    namespace MPL = MaterialPropertyLib;
    MPL::VariableArray variables;

    *out = media_data.medium.property(MPL::PropertyType::biot_coefficient)
               .template value<double>(variables, x_t.x, x_t.t, x_t.dt);
}
}  // namespace ProcessLib::ThermoRichardsMechanics
