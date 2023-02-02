/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "Swelling.h"

#include <Eigen/LU>

#include "MaterialLib/MPL/Utils/FormEigenTensor.h"

namespace ProcessLib::ThermoRichardsMechanics
{
namespace ConstitutiveStress_StrainTemperature
{
template <int DisplacementDim>
void SwellingModel<DisplacementDim>::eval(
    SpaceTimeData const& x_t, MediaData const& media_data,
    ElasticTangentStiffnessData<DisplacementDim> const& C_el_data,
    StrainData<DisplacementDim> const& eps_data,
    PrevState<StrainData<DisplacementDim>> const& eps_prev_data,
    SaturationData const& S_L_data, SaturationDataDeriv const& dS_L_data,
    PrevState<SaturationData> const& S_L_prev_data,
    PrevState<SwellingDataStateful<DisplacementDim>> const& prev_state,
    SwellingDataStateful<DisplacementDim>& state,
    SwellingDataStateless<DisplacementDim>& out) const
{
    namespace MPL = MaterialPropertyLib;

    static constexpr int kelvin_vector_size =
        MathLib::KelvinVector::kelvin_vector_dimensions(DisplacementDim);
    using Invariants = MathLib::KelvinVector::Invariants<kelvin_vector_size>;

    auto const& solid_phase = media_data.solid;

    if (!solid_phase.hasProperty(MPL::PropertyType::swelling_stress_rate))
    {
        out.eps_m.setZero();
        out.J_up_BT_K_N.setZero();
        return;
    }

    MPL::VariableArray variables;
    variables.liquid_saturation = S_L_data.S_L;

    MPL::VariableArray variables_prev;
    variables_prev.liquid_saturation = S_L_prev_data->S_L;

    auto const& identity2 = MathLib::KelvinVector::Invariants<
        MathLib::KelvinVector::kelvin_vector_dimensions(
            DisplacementDim)>::identity2;

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

    auto const C_el_inv = C_el_data.C_el.inverse().eval();

    // !!! Misusing volumetric strain for mechanical volumetric
    // strain just to update the transport porosity !!!
    variables.volumetric_strain =
        Invariants::trace(eps_data.eps) +
        identity2.transpose() * C_el_inv * state.sigma_sw;
    variables_prev.volumetric_strain =
        Invariants::trace(eps_prev_data->eps) +
        identity2.transpose() * C_el_inv * prev_state->sigma_sw;

    out.eps_m.noalias() = C_el_inv * (state.sigma_sw - prev_state->sigma_sw);

    using DimMatrix = Eigen::Matrix<double, 3, 3>;
    auto const dsigma_sw_dS_L =
        MathLib::KelvinVector::tensorToKelvin<DisplacementDim>(
            solid_phase.property(MPL::PropertyType::swelling_stress_rate)
                .template dValue<DimMatrix>(variables, variables_prev,
                                            MPL::Variable::liquid_saturation,
                                            x_t.x, x_t.t, x_t.dt));

    out.J_up_BT_K_N = -dsigma_sw_dS_L * dS_L_data.dS_L_dp_cap;
}

template struct SwellingModel<2>;
template struct SwellingModel<3>;
}  // namespace ConstitutiveStress_StrainTemperature
}  // namespace ProcessLib::ThermoRichardsMechanics
