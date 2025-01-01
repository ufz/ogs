/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include "Base.h"
#include "FluidDensity.h"
#include "Porosity.h"
#include "Saturation.h"
#include "SolidDensity.h"

namespace ProcessLib::TH2M
{
namespace ConstitutiveRelations
{
template <int DisplacementDim>
using VolumetricBodyForce =
    BaseLib::StrongType<GlobalDimVector<DisplacementDim>,
                        struct VolumetricBodyForceTag>;

template <int DisplacementDim>
struct GravityModel
{
    void eval(
        FluidDensityData const& fluid_density_data,
        PorosityData const& porosity_data,
        SaturationData const& S_L_data,
        SolidDensityData const& solid_density_data,
        SpecificBodyForceData<DisplacementDim> const& specific_body_force,
        VolumetricBodyForce<DisplacementDim>& volumetric_body_force) const;
};

extern template struct GravityModel<2>;
extern template struct GravityModel<3>;
}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::TH2M
