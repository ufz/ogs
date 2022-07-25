/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "SolidCompressibility.h"

namespace ProcessLib::ThermoRichardsMechanics
{
template <int DisplacementDim>
void SolidCompressibilityModel<DisplacementDim>::eval(
    const SpaceTimeData& x_t,
    const BiotData& biot_data,
    const ElasticTangentStiffnessData<DisplacementDim>& C_el_data,
    SolidCompressibilityData& out) const
{
    out.beta_SR = (1 - biot_data.alpha) /
                  solid_material_.getBulkModulus(x_t.t, x_t.x, &C_el_data.C_el);
}

template struct SolidCompressibilityModel<2>;
template struct SolidCompressibilityModel<3>;

}  // namespace ProcessLib::ThermoRichardsMechanics
