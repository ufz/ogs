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
namespace ConstitutiveStress_StrainTemperature
{
template <int DisplacementDim>
void SolidMechanicsModel<DisplacementDim>::eval(
    const SpaceTimeData& x_t,
    SolidThermalExpansionData<DisplacementDim> const& s_therm_exp_data,
    SwellingDataStateless<DisplacementDim> const& swelling_data,
    TemperatureData<DisplacementDim> const& T_data,
    CapillaryPressureData<DisplacementDim> const& p_cap_data,
    BiotData const& biot_data,
    BishopsData const& bishops_data,
    SaturationDataDeriv const& dS_L_data,
    StrainData<DisplacementDim> const& eps_data,
    PrevState<StrainData<DisplacementDim>> const& eps_prev_data,
    MaterialStateData<DisplacementDim>& mat_state,
    PrevState<EffectiveStressData<DisplacementDim>> const& sigma_eff_prev_data,
    EffectiveStressData<DisplacementDim>& sigma_eff_data,
    PrevState<MechanicalStrainData<DisplacementDim>> const& eps_m_prev_data,
    MechanicalStrainData<DisplacementDim>& eps_m_data,
    TotalStressData<DisplacementDim>& total_stress_data,
    EquivalentPlasticStrainData& equiv_plast_strain_data,
    SolidMechanicsDataStateless<DisplacementDim>& out) const
{
    namespace MPL = MaterialPropertyLib;
    MPL::VariableArray variables;

    double const T_prev = T_data.T_prev;
    double const dT = T_data.T - T_prev;

    eps_m_data.eps_m.noalias() =
        eps_m_prev_data->eps_m + eps_data.eps - eps_prev_data->eps -
        s_therm_exp_data.solid_linear_thermal_expansivity_vector * dT +
        swelling_data.eps_m;

    variables.mechanical_strain.emplace<KelvinVector<DisplacementDim>>(
        eps_m_data.eps_m);
    variables.temperature = T_data.T;

    MPL::VariableArray variables_prev;
    variables_prev.stress.emplace<KelvinVector<DisplacementDim>>(
        sigma_eff_prev_data->sigma_eff);
    variables_prev.mechanical_strain.emplace<KelvinVector<DisplacementDim>>(
        eps_m_prev_data->eps_m);
    variables_prev.temperature = T_prev;

    auto solution = solid_material_.integrateStress(
        variables_prev, variables, x_t.t, x_t.x, x_t.dt,
        *mat_state.material_state_variables);

    if (!solution)
    {
        OGS_FATAL("Computation of local constitutive relation failed.");
    }

    std::tie(sigma_eff_data.sigma_eff, mat_state.material_state_variables,
             out.stiffness_tensor) = std::move(*solution);

    auto const& identity2 = MathLib::KelvinVector::Invariants<
        MathLib::KelvinVector::kelvin_vector_dimensions(
            DisplacementDim)>::identity2;
    total_stress_data.sigma_total.noalias() =
        sigma_eff_data.sigma_eff +
        biot_data() * bishops_data.chi_S_L * p_cap_data.p_cap * identity2;

    out.J_uT_BT_K_N.noalias() =  // TODO is this thermal stress?
        -out.stiffness_tensor *
        s_therm_exp_data.solid_linear_thermal_expansivity_vector;

    double const J_up_X_BTI2N =
        -biot_data() *
        (bishops_data.chi_S_L +
         bishops_data.dchi_dS_L * p_cap_data.p_cap * dS_L_data.dS_L_dp_cap);

    out.J_up_BT_K_N.noalias() =
        swelling_data.J_up_BT_K_N + J_up_X_BTI2N * identity2;

    equiv_plast_strain_data.equivalent_plastic_strain =
        mat_state.material_state_variables->getEquivalentPlasticStrain();
}

template struct SolidMechanicsModel<2>;
template struct SolidMechanicsModel<3>;

}  // namespace ConstitutiveStress_StrainTemperature
}  // namespace ProcessLib::ThermoRichardsMechanics
