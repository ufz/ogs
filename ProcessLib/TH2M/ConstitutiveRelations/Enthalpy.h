// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "Base.h"
#include "Enthalpy.h"
#include "FluidDensity.h"
#include "PhaseTransitionData.h"
#include "Porosity.h"
#include "ProcessLib/Reflection/ReflectionData.h"
#include "Saturation.h"
#include "SolidDensity.h"
#include "SolidHeatCapacity.h"

namespace ProcessLib::TH2M
{
namespace ConstitutiveRelations
{
struct FluidEnthalpyData
{
    double h_G = nan;
    double h_L = nan;

    static auto reflect()
    {
        using Self = FluidEnthalpyData;
        namespace R = ProcessLib::Reflection;

        return std::tuple{R::makeReflectionData("enthalpy_gas", &Self::h_G),
                          R::makeReflectionData("enthalpy_liquid", &Self::h_L)};
    }
};

struct SolidEnthalpyData
{
    double h_S = nan;

    static auto reflect()
    {
        using Self = SolidEnthalpyData;
        namespace R = ProcessLib::Reflection;

        return std::tuple{R::makeReflectionData("enthalpy_solid", &Self::h_S)};
    }
};

struct EffectiveVolumetricEnthalpy
{
    double rho_h_eff = nan;
};

struct EffectiveVolumetricEnthalpyDerivatives
{
    double drho_h_eff_dT = nan;
    double drho_h_eff_dp_GR = nan;
    double drho_h_eff_dp_cap = nan;
};

struct EffectiveVolumetricEnthalpyModel
{
    void eval(
        FluidDensityData const& fluid_density_data,
        FluidEnthalpyData const& fluid_enthalpy_data,
        PorosityData const& porosity_data,
        SaturationData const& S_L_data,
        SolidDensityData const& solid_density_data,
        SolidEnthalpyData const& solid_enthalpy_data,
        EffectiveVolumetricEnthalpy& effective_volumetric_enthalpy_data) const;

    void dEval(FluidDensityData const& fluid_density_data,
               FluidEnthalpyData const& fluid_enthalpy_data,
               PhaseTransitionData const& phase_transition_data,
               PorosityData const& porosity_data,
               PorosityDerivativeData const& porosity_d_data,
               SaturationData const& S_L_data,
               SolidDensityData const& solid_density_data,
               SolidDensityDerivativeData const& solid_density_d_data,
               SolidEnthalpyData const& solid_enthalpy_data,
               SolidHeatCapacityData const& solid_heat_capacity_data,
               EffectiveVolumetricEnthalpyDerivatives&
                   effective_volumetric_enthalpy_d_data) const;
};

struct SolidEnthalpyModel
{
    void eval(SolidHeatCapacityData const& solid_heat_capacity_data,
              TemperatureData const& T_data,
              SolidEnthalpyData& solid_enthalpy_data) const;
};
}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::TH2M
