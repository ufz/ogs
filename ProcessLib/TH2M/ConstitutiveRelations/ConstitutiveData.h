/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include "Advection.h"
#include "Biot.h"
#include "Bishops.h"
#include "CEquation.h"
#include "ConstitutiveDensity.h"
#include "DarcyVelocity.h"
#include "DiffusionVelocity.h"
#include "ElasticTangentStiffnessData.h"
#include "Enthalpy.h"
#include "EquivalentPlasticStrainData.h"
#include "FluidDensity.h"
#include "InternalEnergy.h"
#include "MassMoleFractions.h"
#include "MechanicalStrain.h"
#include "PermeabilityData.h"
#include "PhaseTransitionData.h"
#include "Porosity.h"
#include "ProcessLib/ConstitutiveRelations/StrainData.h"
#include "ProcessLib/ConstitutiveRelations/StressData.h"
#include "ProcessLib/Reflection/ReflectionData.h"
#include "PureLiquidDensity.h"
#include "Saturation.h"
#include "SolidCompressibility.h"
#include "SolidDensity.h"
#include "SolidHeatCapacity.h"
#include "SolidMechanics.h"
#include "SolidThermalExpansion.h"
#include "Swelling.h"
#include "ThermalConductivity.h"
#include "TotalStress.h"
#include "VapourPartialPressure.h"
#include "Viscosity.h"
#include "WEquation.h"

namespace ProcessLib::TH2M
{
namespace ConstitutiveRelations
{
/// Data whose state must be tracked by the process.
template <int DisplacementDim>
struct StatefulData
{
    SaturationData S_L_data;
    SwellingDataStateful<DisplacementDim> swelling_data;
    ProcessLib::ConstitutiveRelations::StressData<DisplacementDim>
        eff_stress_data;
    MechanicalStrainData<DisplacementDim> mechanical_strain_data;
    PureLiquidDensityData rho_W_LR;
    ConstituentDensityData constituent_density_data;
    InternalEnergyData internal_energy_data;

    static auto reflect()
    {
        using Self = StatefulData<DisplacementDim>;

        return Reflection::reflectWithoutName(
            &Self::S_L_data, &Self::swelling_data, &Self::eff_stress_data);
    }
};

template <int DisplacementDim>
struct StatefulDataPrev
{
    PrevState<SaturationData> S_L_data;
    PrevState<SwellingDataStateful<DisplacementDim>> swelling_data;
    PrevState<ProcessLib::ConstitutiveRelations::StressData<DisplacementDim>>
        eff_stress_data;
    PrevState<MechanicalStrainData<DisplacementDim>> mechanical_strain_data;
    PrevState<PureLiquidDensityData> rho_W_LR;
    PrevState<ConstituentDensityData> constituent_density_data;
    PrevState<InternalEnergyData> internal_energy_data;

    StatefulDataPrev<DisplacementDim>& operator=(
        StatefulData<DisplacementDim> const& state)
    {
        S_L_data = state.S_L_data;
        swelling_data = state.swelling_data;
        eff_stress_data = state.eff_stress_data;
        mechanical_strain_data = state.mechanical_strain_data;
        rho_W_LR = state.rho_W_LR;
        constituent_density_data = state.constituent_density_data;
        internal_energy_data = state.internal_energy_data;

        return *this;
    }
};

/// Data that is needed for output purposes.
template <int DisplacementDim>
struct OutputData
{
    ProcessLib::ConstitutiveRelations::StrainData<DisplacementDim> eps_data;
    PermeabilityData<DisplacementDim> permeability_data;
    EnthalpyData enthalpy_data;
    MassMoleFractionsData mass_mole_fractions_data;
    FluidDensityData fluid_density_data;
    VapourPartialPressureData vapour_pressure_data;
    PorosityData porosity_data;
    SolidDensityData solid_density_data;
    DiffusionVelocityData<DisplacementDim> diffusion_velocity_data;
    DarcyVelocityData<DisplacementDim> darcy_velocity_data;

    static auto reflect()
    {
        using Self = OutputData<DisplacementDim>;

        return Reflection::reflectWithoutName(&Self::eps_data,
                                              &Self::permeability_data,
                                              &Self::enthalpy_data,
                                              &Self::mass_mole_fractions_data,
                                              &Self::fluid_density_data,
                                              &Self::vapour_pressure_data,
                                              &Self::porosity_data,
                                              &Self::solid_density_data,
                                              &Self::diffusion_velocity_data,
                                              &Self::darcy_velocity_data);
    }
};

/// Data that is needed for the equation system assembly.
template <int DisplacementDim>
struct ConstitutiveData
{
    SolidMechanicsDataStateless<DisplacementDim> s_mech_data;
};

/// Data that stores intermediate values, which are not needed outside the
/// constitutive setting.
template <int DisplacementDim>
struct ConstitutiveTempData
{
    SwellingDataStateless<DisplacementDim> swelling_data;
    ElasticTangentStiffnessData<DisplacementDim> C_el_data;
    BiotData biot_data;
    SolidCompressibilityData beta_p_SR;
    SaturationDataDeriv dS_L_dp_cap;
    BishopsData chi_S_L;
    SolidThermalExpansionData<DisplacementDim> s_therm_exp_data;
    TotalStressData<DisplacementDim> total_stress_data;
    EquivalentPlasticStrainData equivalent_plastic_strain_data;
    ViscosityData viscosity_data;
    PhaseTransitionData phase_transition_data;
    PorosityDerivativeData porosity_d_data;
    SolidHeatCapacityData solid_heat_capacity_data;
    ThermalConductivityData<DisplacementDim> thermal_conductivity_data;
    EffectiveVolumetricEnthalpy effective_volumetric_enthalpy_data;
    AdvectionData<DisplacementDim> advection_data;
    FC1Data<DisplacementDim> fC_1;
    FC2aData fC_2a;
    FC3aData fC_3a;
    FC4LCpGData<DisplacementDim> fC_4_LCpG;
    FC4LCpCData<DisplacementDim> fC_4_LCpC;
    FC4LCTData<DisplacementDim> fC_4_LCT;
    FC4MCpGData fC_4_MCpG;
    FC4MCpCData fC_4_MCpC;
    FC4MCTData fC_4_MCT;
    FC4MCuData fC_4_MCu;

    FW1Data<DisplacementDim> fW_1;

    using DisplacementDimVector = Eigen::Matrix<double, DisplacementDim, 1>;
    using DisplacementDimMatrix =
        Eigen::Matrix<double, DisplacementDim, DisplacementDim>;

    DisplacementDimVector drho_GR_h_w_eff_dp_GR_Npart;
    DisplacementDimMatrix drho_GR_h_w_eff_dp_GR_gradNpart;
    DisplacementDimVector drho_LR_h_w_eff_dp_cap_Npart;
    DisplacementDimMatrix drho_LR_h_w_eff_dp_cap_gradNpart;
    DisplacementDimVector drho_GR_h_w_eff_dT;
    DisplacementDimMatrix dfW_4_LWpG_a_dp_GR;
    DisplacementDimMatrix dfW_4_LWpG_a_dp_cap;
    DisplacementDimMatrix dfW_4_LWpG_a_dT;
    DisplacementDimMatrix dfW_4_LWpG_d_dp_GR;
    DisplacementDimMatrix dfW_4_LWpG_d_dp_cap;
    DisplacementDimMatrix dfW_4_LWpG_d_dT;
    DisplacementDimMatrix dfW_4_LWpC_a_dp_GR;
    DisplacementDimMatrix dfW_4_LWpC_a_dp_cap;
    DisplacementDimMatrix dfW_4_LWpC_a_dT;
    DisplacementDimMatrix dfW_4_LWpC_d_dp_GR;
    DisplacementDimMatrix dfW_4_LWpC_d_dp_cap;
    DisplacementDimMatrix dfW_4_LWpC_d_dT;
    DisplacementDimMatrix dk_over_mu_G_dp_cap;
    DisplacementDimMatrix dk_over_mu_L_dp_cap;
    double dfW_2a_dp_GR = std::numeric_limits<double>::quiet_NaN();
    double dfW_2b_dp_GR = std::numeric_limits<double>::quiet_NaN();
    double dfW_2a_dp_cap = std::numeric_limits<double>::quiet_NaN();
    double dfW_2b_dp_cap = std::numeric_limits<double>::quiet_NaN();
    double dfW_2a_dT = std::numeric_limits<double>::quiet_NaN();
    double dfW_2b_dT = std::numeric_limits<double>::quiet_NaN();
    double dfW_3a_dp_GR = std::numeric_limits<double>::quiet_NaN();
    double dfW_3a_dp_cap = std::numeric_limits<double>::quiet_NaN();
    double dfW_3a_dT = std::numeric_limits<double>::quiet_NaN();
};

/// Data that stores intermediate values (derivatives), which are not needed
/// outside the constitutive setting.
template <int DisplacementDim>
struct DerivativesData
{
    AdvectionDerivativeData<DisplacementDim> advection_d_data;
    ThermalConductivityDerivativeData<DisplacementDim>
        thermal_conductivity_d_data;
    SolidDensityDerivativeData solid_density_d_data;
    EffectiveVolumetricInternalEnergyDerivatives
        effective_volumetric_internal_energy_d_data;
    EffectiveVolumetricEnthalpyDerivatives effective_volumetric_enthalpy_d_data;
    FC2aDerivativeData dfC_2a;
    FC3aDerivativeData dfC_3a;
    FC4LCpGDerivativeData<DisplacementDim> dfC_4_LCpG;
    FC4LCpCDerivativeData<DisplacementDim> dfC_4_LCpC;
    FC4MCpGDerivativeData dfC_4_MCpG;
    FC4MCTDerivativeData dfC_4_MCT;
    FC4MCuDerivativeData dfC_4_MCu;
};

}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::TH2M
