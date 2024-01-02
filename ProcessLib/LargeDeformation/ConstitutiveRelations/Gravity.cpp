/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "Gravity.h"

namespace ProcessLib::LargeDeformation
{
template <int DisplacementDim>
void GravityModel<DisplacementDim>::eval(
    SolidDensity const& rho_SR, VolumetricBodyForce<DisplacementDim>& out) const
{
    *out = *rho_SR * specific_body_force_;
}

template struct GravityModel<2>;
template struct GravityModel<3>;
}  // namespace ProcessLib::LargeDeformation
