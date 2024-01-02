/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "SolidMechanics.h"

namespace ProcessLib::LargeDeformation
{
namespace ConstitutiveRelations
{
template <int DisplacementDim>
void SolidMechanicsModel<DisplacementDim>::eval(
    SpaceTimeData const& x_t,
    Temperature const& temperature,
    DeformationGradientData<DisplacementDim> const& deformation_gradient_data,
    PrevState<DeformationGradientData<DisplacementDim>> const&
        deformation_gradient_data_prev,
    MaterialStateData<DisplacementDim>& mat_state,
    PrevState<StressData<DisplacementDim>> const& stress_data_prev,
    StressData<DisplacementDim>& stress_data,
    SolidMechanicsDataStateless<DisplacementDim>& current_stateless) const
{
    namespace MPL = MaterialPropertyLib;

    // current state
    MPL::VariableArray variables;
    {
        // thermodynamic forces
        variables.stress = stress_data.sigma;

        // gradient
        variables.deformation_gradient =
            deformation_gradient_data.deformation_gradient;

        // external state variables
        variables.temperature = *temperature;
    }

    // previous state
    MPL::VariableArray variables_prev;
    {
        // thermodynamic forces
        variables_prev.stress = stress_data_prev->sigma;

        // gradient
        variables_prev.deformation_gradient =
            deformation_gradient_data_prev->deformation_gradient;

        // external state variables
        variables_prev.temperature = *temperature;
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

    stress_data.sigma =
        view.block(MSM::second_piola_kirchhoff_stress, tdyn_forces_data);
    mat_state.material_state_variables = std::move(std::get<1>(*solution));

    auto const& tangent_operator_data = std::get<2>(*solution);

    current_stateless.stiffness_tensor = tangent_operator_blocks_view_.block(
        MSM::second_piola_kirchhoff_stress, MSM::green_lagrange_strain,
        tangent_operator_data);
}

template struct SolidMechanicsModel<2>;
template struct SolidMechanicsModel<3>;
}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::LargeDeformation
