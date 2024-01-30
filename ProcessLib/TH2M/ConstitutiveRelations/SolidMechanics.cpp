/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "SolidMechanics.h"

namespace ProcessLib::TH2M
{
namespace ConstitutiveRelations
{
template <int DisplacementDim>
void SolidMechanicsModel<DisplacementDim>::eval(
    const SpaceTimeData& x_t,
    TemperatureData const& T_data,
    MechanicalStrainData<DisplacementDim> const& mechanical_strain_data,
    PrevState<MechanicalStrainData<DisplacementDim>> const&
        mechanical_strain_prev_data,
    PrevState<
        ProcessLib::ConstitutiveRelations::StressData<DisplacementDim>> const&
        eff_stress_prev_data,
    ProcessLib::ConstitutiveRelations::StressData<DisplacementDim>&
        eff_stress_data,
    MaterialStateData<DisplacementDim>& mat_state,
    SolidMechanicsDataStateless<DisplacementDim>& out,
    EquivalentPlasticStrainData& equivalent_plastic_strain) const
{
    namespace MPL = MaterialPropertyLib;
    MPL::VariableArray variables;

    variables.mechanical_strain.emplace<KelvinVector<DisplacementDim>>(
        mechanical_strain_data.eps_m);
    variables.temperature = T_data.T;

    MPL::VariableArray variables_prev;
    variables_prev.stress.emplace<KelvinVector<DisplacementDim>>(
        eff_stress_prev_data->sigma);
    variables_prev.mechanical_strain.emplace<KelvinVector<DisplacementDim>>(
        mechanical_strain_prev_data->eps_m);
    variables_prev.temperature = T_data.T_prev;

    auto solution = solid_material_.integrateStress(
        variables_prev, variables, x_t.t, x_t.x, x_t.dt,
        *mat_state.material_state_variables);

    if (!solution)
    {
        OGS_FATAL("Computation of local constitutive relation failed.");
    }

    std::tie(eff_stress_data.sigma, mat_state.material_state_variables,
             out.stiffness_tensor) = std::move(*solution);

    *equivalent_plastic_strain =
        mat_state.material_state_variables->getEquivalentPlasticStrain();
}

template struct SolidMechanicsModel<2>;
template struct SolidMechanicsModel<3>;
}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::TH2M
