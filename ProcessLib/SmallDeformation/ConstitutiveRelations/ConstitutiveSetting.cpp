/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "ConstitutiveSetting.h"

#include "Invoke.h"

namespace ProcessLib::SmallDeformation
{
namespace ConstitutiveRelations
{
template <int DisplacementDim>
void ConstitutiveSetting<DisplacementDim>::eval(
    ConstitutiveModels<DisplacementDim>& models, double const t,
    double const dt, ParameterLib::SpatialPosition const& x_position,
    MaterialPropertyLib::Medium const& medium, double const T_ref,
    KelvinVector<DisplacementDim> const& eps,
    KelvinVector<DisplacementDim> const& eps_prev,
    StatefulData<DisplacementDim>& state,
    StatefulDataPrev<DisplacementDim> const& prev_state,
    MaterialStateData<DisplacementDim>& mat_state,
    ConstitutiveTempData<DisplacementDim>& tmp,
    OutputData<DisplacementDim>& out,
    ConstitutiveData<DisplacementDim>& cd) const
{
    namespace MPL = MaterialPropertyLib;

    auto& eps_data = out.eps_data;
    eps_data.eps = eps;
    auto& eps_data_prev = tmp.eps_data_prev;
    eps_data_prev->eps = eps_prev;
    auto& rho_SR = tmp.rho_SR;

    auto& s_mech_data = cd.s_mech_data;
    auto& volumetric_body_force = cd.volumetric_body_force;

    auto& free_energy_density_data = out.free_energy_density_data;

    Temperature const T{T_ref};
    SpaceTimeData const x_t{x_position, t, dt};
    MediaData const media_data{medium};

    assertEvalArgsUnique(models.s_mech_model);
    models.s_mech_model.eval(x_t, T, eps_data, eps_data_prev, mat_state,
                             prev_state.stress_data, state.stress_data,
                             s_mech_data, free_energy_density_data);

    assertEvalArgsUnique(models.rho_S_model);
    models.rho_S_model.eval(x_t, media_data, T, rho_SR);

    assertEvalArgsUnique(models.gravity_model);
    models.gravity_model.eval(rho_SR, volumetric_body_force);
}

template struct ConstitutiveSetting<2>;
template struct ConstitutiveSetting<3>;
}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::SmallDeformation
