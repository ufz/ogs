// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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
