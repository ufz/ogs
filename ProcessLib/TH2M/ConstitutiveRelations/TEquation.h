// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "Base.h"
#include "ConstitutiveDensity.h"
#include "DarcyVelocity.h"
#include "DiffusionVelocity.h"
#include "Enthalpy.h"
#include "FluidDensity.h"
#include "InternalEnergy.h"
#include "PermeabilityData.h"
#include "Viscosity.h"

namespace ProcessLib::TH2M
{
namespace ConstitutiveRelations
{
struct FT1Data
{
    double m = nan;
};

struct FT1DerivativeData
{
    double dp_GR = nan;
    double dp_cap = nan;
    double dT = nan;
};

struct FT1Model
{
    void eval(double const dt,
              InternalEnergyData const& internal_energy_data,
              PrevState<InternalEnergyData> const& internal_energy_data_prev,
              FT1Data& fT_1) const;

    void dEval(double const dt,
               EffectiveVolumetricInternalEnergyDerivatives const&
                   effective_volumetric_internal_energy_d_data,
               FT1DerivativeData& dfT_1) const;
};

template <int DisplacementDim>
struct FT2Data
{
    GlobalDimVector<DisplacementDim> A;
};

template <int DisplacementDim>
struct FT2DerivativeData
{
    GlobalDimVector<DisplacementDim> dp_GR_Npart;
    GlobalDimMatrix<DisplacementDim> dp_GR_gradNpart;
    GlobalDimVector<DisplacementDim> dp_cap_Npart;
    GlobalDimMatrix<DisplacementDim> dp_cap_gradNpart;
    GlobalDimVector<DisplacementDim> dT;
};

template <int DisplacementDim>
struct FT2Model
{
    void eval(DarcyVelocityData<DisplacementDim> const& darcy_velocity_data,
              FluidDensityData const& fluid_density_data,
              FluidEnthalpyData const& fluid_enthalpy_data,
              FT2Data<DisplacementDim>& fT_2) const;

    void dEval(
        DarcyVelocityData<DisplacementDim> const& darcy_velocity_data,
        FluidDensityData const& fluid_density_data,
        FluidEnthalpyData const& fluid_enthalpy_data,
        PermeabilityData<DisplacementDim> const& permeability_data,
        PhaseTransitionData const& phase_transition_data,
        SpecificBodyForceData<DisplacementDim> const& specific_body_force,
        ViscosityData const& viscosity_data,
        FT2DerivativeData<DisplacementDim>& dfT_2) const;
};

extern template struct FT2Model<2>;
extern template struct FT2Model<3>;

template <int DisplacementDim>
struct FT3Data
{
    double N = nan;
    GlobalDimVector<DisplacementDim> gradN;
};

template <int DisplacementDim>
struct FT3Model
{
    void eval(
        ConstituentDensityData const& constituent_density_data,
        DarcyVelocityData<DisplacementDim> const& darcy_velocity_data,
        DiffusionVelocityData<DisplacementDim> const& diffusion_velocity_data,
        FluidDensityData const& fluid_density_data,
        PhaseTransitionData const& phase_transition_data,
        SpecificBodyForceData<DisplacementDim> const& specific_body_force,
        FT3Data<DisplacementDim>& fT_3) const;
};

extern template struct FT3Model<2>;
extern template struct FT3Model<3>;
}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::TH2M
