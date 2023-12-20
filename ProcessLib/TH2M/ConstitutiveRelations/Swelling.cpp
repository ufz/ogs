/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "Swelling.h"

#include <Eigen/LU>

#include "MaterialLib/MPL/Utils/FormEigenTensor.h"

namespace ProcessLib::TH2M
{
namespace ConstitutiveRelations
{
template <int DisplacementDim>
void SwellingModel<DisplacementDim>::eval(
    SpaceTimeData const& x_t, MediaData const& media_data,
    ElasticTangentStiffnessData<DisplacementDim> const& C_el_data,
    SaturationData const& S_L_data,
    PrevState<SaturationData> const& S_L_prev_data,
    PrevState<SwellingDataStateful<DisplacementDim>> const& prev_state,
    SwellingDataStateful<DisplacementDim>& state,
    SwellingDataStateless<DisplacementDim>& out) const
{
    namespace MPL = MaterialPropertyLib;

    auto const& solid_phase = media_data.solid;

    if (!solid_phase.hasProperty(MPL::PropertyType::swelling_stress_rate))
    {
        out.eps_m.setZero();
        state.sigma_sw.setZero();
        return;
    }

    MPL::VariableArray variables;
    variables.liquid_saturation = S_L_data.S_L;

    MPL::VariableArray variables_prev;
    variables_prev.liquid_saturation = S_L_prev_data->S_L;

    // Swelling and possibly volumetric strain rate update.

    // If there is swelling, compute it. Update volumetric strain rate,
    // s.t. it corresponds to the mechanical part only.
    state.sigma_sw = prev_state->sigma_sw;

    auto const sigma_sw_dot =
        MathLib::KelvinVector::tensorToKelvin<DisplacementDim>(
            MPL::formEigenTensor<3>(
                solid_phase[MPL::PropertyType::swelling_stress_rate].value(
                    variables, variables_prev, x_t.x, x_t.t, x_t.dt)));
    state.sigma_sw += sigma_sw_dot * x_t.dt;

    auto const C_el_inv = C_el_data.stiffness_tensor.inverse().eval();

    out.eps_m.noalias() = C_el_inv * (state.sigma_sw - prev_state->sigma_sw);
}

template struct SwellingModel<2>;
template struct SwellingModel<3>;
}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::TH2M
