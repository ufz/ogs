/**
 * \file
 * \copyright
 * Copyright (c) 2012-2023, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "Gravity.h"

namespace ProcessLib::SmallDeformation
{
template <int DisplacementDim>
void GravityModel<DisplacementDim>::eval(
    SolidDensity const& rho_SR, GravityData<DisplacementDim>& out) const
{
    auto const b = specific_body_force_;

    out.volumetric_body_force = *rho_SR * b;
}

template struct GravityModel<2>;
template struct GravityModel<3>;
}  // namespace ProcessLib::SmallDeformation
