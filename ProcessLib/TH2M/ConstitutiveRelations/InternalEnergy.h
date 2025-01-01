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
#include "BaseLib/StrongType.h"
#include "Enthalpy.h"
#include "FluidDensity.h"
#include "PhaseTransitionData.h"
#include "Porosity.h"
#include "Saturation.h"
#include "SolidDensity.h"
#include "SolidHeatCapacity.h"

namespace ProcessLib::TH2M
{
namespace ConstitutiveRelations
{
using InternalEnergyData =
    BaseLib::StrongType<double, struct InternalEnergyTag>;

struct EffectiveVolumetricInternalEnergyDerivatives
{
    double drho_u_eff_dT = nan;
    double drho_u_eff_dp_GR = nan;
    double drho_u_eff_dp_cap = nan;
};

struct InternalEnergyModel
{
    void eval(FluidDensityData const& fluid_density_data,
              PhaseTransitionData const& phase_transition_data,
              PorosityData const& porosity_data,
              SaturationData const& S_L_data,
              SolidDensityData const& solid_density_data,
              SolidEnthalpyData const& solid_enthalpy_data,
              InternalEnergyData& internal_energy_data) const;

    void dEval(FluidDensityData const& fluid_density_data,
               PhaseTransitionData const& phase_transition_data,
               PorosityData const& porosity_data,
               PorosityDerivativeData const& porosity_d_data,
               SaturationData const& S_L_data,
               SolidDensityData const& solid_density_data,
               SolidDensityDerivativeData const& solid_density_d_data,
               SolidEnthalpyData const& solid_enthalpy_data,
               SolidHeatCapacityData const& solid_heat_capacity_data,
               EffectiveVolumetricInternalEnergyDerivatives&
                   effective_volumetric_internal_energy_d_data) const;
};
}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::TH2M
