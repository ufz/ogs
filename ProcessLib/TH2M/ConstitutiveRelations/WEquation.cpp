/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "WEquation.h"

namespace ProcessLib::TH2M
{
namespace ConstitutiveRelations
{
template <int DisplacementDim>
void FW1Model<DisplacementDim>::eval(
    AdvectionData<DisplacementDim> const& advection_data,
    FluidDensityData const& fluid_density_data,
    FW1Data<DisplacementDim>& fW_1) const
{
    fW_1.A = advection_data.advection_W_G * fluid_density_data.rho_GR +
             advection_data.advection_W_L * fluid_density_data.rho_LR;
}

template struct FW1Model<2>;
template struct FW1Model<3>;
}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::TH2M
