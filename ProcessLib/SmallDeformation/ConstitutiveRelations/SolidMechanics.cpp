/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "SolidMechanics.h"

namespace ProcessLib::SmallDeformation
{
namespace ConstitutiveRelations
{
template <int DisplacementDim>
void SolidMechanicsModel<DisplacementDim>::eval(
    SpaceTimeData const& x_t,
    Temperature const& temperature,
    StrainData<DisplacementDim> const& eps_data,
    PrevState<StrainData<DisplacementDim>> const& eps_data_prev,
    MaterialStateData<DisplacementDim>& mat_state,
    PrevState<StressData<DisplacementDim>> const& stress_data_prev,
    StressData<DisplacementDim>& stress_data,
    SolidMechanicsDataStateless<DisplacementDim>& current_stateless,
    FreeEnergyDensityData& free_energy_density_data) const
{
    namespace MPL = MaterialPropertyLib;

    // current state
    MPL::VariableArray variables;
    {
        // thermodynamic forces
        variables.stress = stress_data.sigma;
        variables.mechanical_strain = eps_data.eps;

        // external state variables
        variables.temperature = *temperature;
    }

    // previous state
    MPL::VariableArray variables_prev;
    {
        // thermodynamic forces
        variables_prev.stress = stress_data_prev->sigma;
        variables_prev.mechanical_strain = eps_data_prev->eps;

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

    std::tie(stress_data.sigma, mat_state.material_state_variables,
             current_stateless.stiffness_tensor) = std::move(*solution);

    free_energy_density_data.free_energy_density =
        solid_material_.computeFreeEnergyDensity(
            x_t.t, x_t.x, x_t.dt, eps_data.eps, stress_data.sigma,
            *mat_state.material_state_variables);
}

template struct SolidMechanicsModel<2>;
template struct SolidMechanicsModel<3>;
}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::SmallDeformation
