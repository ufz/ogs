/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "SolidMechanics.h"

namespace ProcessLib::ThermoRichardsMechanics
{
namespace ConstitutiveStressSaturation_StrainPressureTemperature
{
template <int DisplacementDim>
void SolidMechanicsModel<DisplacementDim>::eval(
    SpaceTimeData const& x_t, TemperatureData<DisplacementDim> const& T_data,
    CapillaryPressureData<DisplacementDim> const& p_cap_data,
    StrainData<DisplacementDim> const& eps_data,
    PrevState<StrainData<DisplacementDim>> const& eps_prev_data,
    MaterialStateData<DisplacementDim>& mat_state,
    PrevState<TotalStressData<DisplacementDim>> const& total_stress_data_prev,
    TotalStressData<DisplacementDim>& total_stress_data,
    EquivalentPlasticStrainData& equiv_plast_strain_data,
    SolidMechanicsDataStateless<DisplacementDim>& current_stateless,
    PrevState<SaturationData> const& S_L_prev_data, SaturationData& S_L_data,
    SaturationDataDeriv& dS_L_data) const
{
    namespace MPL = MaterialPropertyLib;

    double const T_prev = T_data.T_prev;
    auto const& eps_total = eps_data.eps;
    auto const& eps_total_prev = eps_prev_data->eps;
    auto const& sigma_total_prev = total_stress_data_prev->sigma_total;

    // current state
    MPL::VariableArray variables;
    {
        // gradients
        // TODO currently we always pass strain via mechanical_strain
        variables.mechanical_strain = eps_total;
        variables.liquid_phase_pressure = -p_cap_data.p_cap;

        // external state variables
        variables.temperature = T_data.T;
    }

    // previous state
    MPL::VariableArray variables_prev;
    {
        // thermodynamic forces
        variables_prev.stress = sigma_total_prev;
        variables_prev.liquid_saturation = S_L_prev_data->S_L;

        // gradients
        // TODO currently we always pass strain via mechanical_strain
        variables_prev.mechanical_strain = eps_total_prev;
        variables_prev.liquid_phase_pressure = -p_cap_data.p_cap_prev;

        // external state variables
        variables_prev.temperature = T_prev;
    }

    auto solution = solid_material_.integrateStress(
        variables_prev, variables, x_t.t, x_t.x, x_t.dt,
        *mat_state.material_state_variables);

    if (!solution)
    {
        OGS_FATAL("Computation of local constitutive relation failed.");
    }

    auto& tdyn_forces_data = std::get<0>(*solution);

    auto const view = solid_material_.createThermodynamicForcesView();

    total_stress_data.sigma_total = view.block(MSM::stress, tdyn_forces_data);
    S_L_data.S_L = view.block(MSM::saturation, tdyn_forces_data);
    mat_state.material_state_variables = std::move(std::get<1>(*solution));

    auto const& tangent_operator_data = std::get<2>(*solution);

    current_stateless.stiffness_tensor = tangent_operator_blocks_view_.block(
        MSM::stress, MSM::strain, tangent_operator_data);

    dS_L_data.dS_L_dp_cap = -tangent_operator_blocks_view_.block(
        MSM::saturation, MSM::liquid_pressure, tangent_operator_data);

    current_stateless.J_uT_BT_K_N = tangent_operator_blocks_view_.block(
        MSM::stress, MSM::temperature, tangent_operator_data);

    current_stateless.J_up_BT_K_N = tangent_operator_blocks_view_.block(
        MSM::stress, MSM::liquid_pressure, tangent_operator_data);

    equiv_plast_strain_data.equivalent_plastic_strain =
        mat_state.material_state_variables->getEquivalentPlasticStrain();
}

template struct SolidMechanicsModel<2>;
template struct SolidMechanicsModel<3>;
}  // namespace ConstitutiveStressSaturation_StrainPressureTemperature
}  // namespace ProcessLib::ThermoRichardsMechanics
