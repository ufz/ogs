/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include "Advection.h"
#include "Biot.h"
#include "Bishops.h"
#include "CEquation.h"
#include "DiffusionVelocity.h"
#include "ElasticTangentStiffnessModel.h"
#include "Enthalpy.h"
#include "Gravity.h"
#include "InternalEnergy.h"
#include "MechanicalStrainModel.h"
#include "PermeabilityModel.h"
#include "PhaseTransitionModel.h"
#include "Porosity.h"
#include "PureLiquidDensity.h"
#include "Saturation.h"
#include "SolidCompressibility.h"
#include "SolidDensity.h"
#include "SolidHeatCapacity.h"
#include "SolidMechanics.h"
#include "SolidThermalExpansion.h"
#include "Swelling.h"
#include "TEquation.h"
#include "ThermalConductivity.h"
#include "TotalStress.h"
#include "TransportPorosity.h"
#include "UEquation.h"
#include "Viscosity.h"
#include "WEquation.h"

namespace ProcessLib::TH2M
{
namespace ConstitutiveRelations
{
/// Constitutive models used for assembly.
template <int DisplacementDim>
struct ConstitutiveModels
{
    explicit ConstitutiveModels(
        SolidConstitutiveRelation<DisplacementDim> const& solid_material,
        PhaseTransitionModel const& phase_transition_model)
        : elastic_tangent_stiffness_model(solid_material),
          beta_p_SR_model(solid_material),
          s_mech_model(solid_material),
          phase_transition_model(phase_transition_model)
    {
    }

    ElasticTangentStiffnessModel<DisplacementDim>
        elastic_tangent_stiffness_model;
    BiotModel biot_model;
    SolidCompressibilityModel<DisplacementDim,
                              SolidConstitutiveRelation<DisplacementDim>>
        beta_p_SR_model;
    SaturationModel S_L_model;
    BishopsModel chi_S_L_model;
    SwellingModel<DisplacementDim> swelling_model;
    SolidThermalExpansionModel<DisplacementDim> s_therm_exp_model;
    MechanicalStrainModel<DisplacementDim> mechanical_strain_model;
    SolidMechanicsModel<DisplacementDim> s_mech_model;
    TotalStressModel<DisplacementDim> total_stress_model;
    PermeabilityModel<DisplacementDim> permeability_model;
    PureLiquidDensityModel pure_liquid_density_model;
    PhaseTransitionModel const& phase_transition_model;
    ViscosityModel viscosity_model;
    BishopsPrevModel chi_S_L_prev_model;
    PorosityModel<DisplacementDim> porosity_model;
    TransportPorosityModel<DisplacementDim> transport_porosity_model;
    SolidDensityModel<DisplacementDim> solid_density_model;
    SolidHeatCapacityModel solid_heat_capacity_model;
    ThermalConductivityModel<DisplacementDim> thermal_conductivity_model;
    AdvectionModel<DisplacementDim> advection_model;
    InternalEnergyModel internal_energy_model;
    EffectiveVolumetricEnthalpyModel effective_volumetric_enthalpy_model;
    GravityModel<DisplacementDim> gravity_model;
    DiffusionVelocityModel<DisplacementDim> diffusion_velocity_model;
    SolidEnthalpyModel solid_enthalpy_model;
    DarcyVelocityModel<DisplacementDim> darcy_velocity_model;
    FC1Model<DisplacementDim> fC_1_model;
    FC2aModel fC_2a_model;
    FC3aModel fC_3a_model;
    FC4LCpGModel<DisplacementDim> fC_4_LCpG_model;
    FC4LCpCModel<DisplacementDim> fC_4_LCpC_model;
    FC4LCTModel<DisplacementDim> fC_4_LCT_model;
    FC4MCpGModel fC_4_MCpG_model;
    FC4MCpCModel fC_4_MCpC_model;
    FC4MCTModel<DisplacementDim> fC_4_MCT_model;
    FC4MCuModel fC_4_MCu_model;

    FW1Model<DisplacementDim> fW_1_model;
    FW2Model fW_2_model;
    FW3aModel fW_3a_model;
    FW4LWpGModel<DisplacementDim> fW_4_LWpG_model;
    FW4LWpCModel<DisplacementDim> fW_4_LWpC_model;
    FW4LWTModel<DisplacementDim> fW_4_LWT_model;
    FW4MWpGModel fW_4_MWpG_model;
    FW4MWpCModel fW_4_MWpC_model;
    FW4MWTModel<DisplacementDim> fW_4_MWT_model;
    FW4MWuModel fW_4_MWu_model;

    FT1Model fT_1_model;
    FT2Model<DisplacementDim> fT_2_model;
    FT3Model<DisplacementDim> fT_3_model;

    FU1KUTModel<DisplacementDim> fu_1_KuT_model;
    FU2KUpCModel fu_2_KupC_model;
};
}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::TH2M
