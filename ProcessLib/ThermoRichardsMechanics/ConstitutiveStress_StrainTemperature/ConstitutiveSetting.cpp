/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ConstitutiveSetting.h"

#include "ProcessLib/Graph/Apply.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/Invoke.h"

namespace ProcessLib::ThermoRichardsMechanics
{
namespace ConstitutiveStress_StrainTemperature
{
template <int DisplacementDim>
void ConstitutiveSetting<DisplacementDim>::eval(
    ConstitutiveModels<DisplacementDim>& models, double const t,
    double const dt, ParameterLib::SpatialPosition const& x_position,
    MaterialPropertyLib::Medium const& medium,
    TemperatureData<DisplacementDim> const& T_data,
    CapillaryPressureData<DisplacementDim> const& p_cap_data,
    KelvinVector<DisplacementDim> const& eps_arg,
    StatefulData<DisplacementDim>& state,
    StatefulDataPrev<DisplacementDim> const& prev_state,
    MaterialStateData<DisplacementDim>& mat_state,
    ConstitutiveTempData<DisplacementDim>& tmp,
    OutputData<DisplacementDim>& out,
    ConstitutiveData<DisplacementDim>& cd) const
{
    namespace MPL = MaterialPropertyLib;

    auto& biot_data = std::get<BiotData>(tmp);
    auto& bishops_data_prev = std::get<PrevState<BishopsData>>(tmp);
    auto& poro_data = std::get<PorosityData>(state);

    SpaceTimeData const x_t{x_position, t, dt};
    MediaData const media_data{medium};

    auto const aux_data = std::tuple{SpaceTimeData{x_position, t, dt},
                                     MediaData{medium}, T_data, p_cap_data};

    auto const mat_state_tuple = std::tie(mat_state);

    namespace G = ProcessLib::Graph;

    // TODO will eps lag one iteration behind? (since it's not updated after
    // solving the global equation system)
    std::get<StrainData<DisplacementDim>>(state).eps.noalias() = eps_arg;

    G::apply(models.elastic_tangent_stiffness_model, aux_data, tmp);
    G::apply(models.biot_model, aux_data, tmp);
    G::apply(models.solid_compressibility_model, aux_data, tmp);
    G::apply(models.S_L_model, aux_data, state, tmp);

    G::apply(models.bishops_model, aux_data, state, tmp);

    assertEvalArgsUnique(models.bishops_model);
    // TODO why not ordinary state tracking?
    models.bishops_model.eval(x_t, media_data,
                              *std::get<PrevState<SaturationData>>(prev_state),
                              *bishops_data_prev);

    G::apply(models.poro_model, aux_data, tmp, state, prev_state);

    if (biot_data() < poro_data.phi)
    {
        OGS_FATAL(
            "ThermoRichardsMechanics: Biot-coefficient {} is smaller than "
            "porosity {} in element/integration point {}/{}.",
            biot_data(), poro_data.phi, *x_position.getElementID(),
            *x_position.getIntegrationPoint());
    }

    G::apply(models.swelling_model, aux_data, state, prev_state, tmp);
    G::apply(models.s_therm_exp_model, aux_data, tmp);
    G::apply(models.s_mech_model, aux_data, tmp, state, prev_state,
             mat_state_tuple, cd);

    G::apply(models.rho_L_model, aux_data, out);

    /* {
        double const p_FR = -bishops_data.chi_S_L * p_cap_data.p_cap;
        // p_SR
        // TODO used by no MPL model
        variables.solid_grain_pressure =
            p_FR -
            Invariants::trace(std::get<s_mech_data>(state).sigma_eff) / (3 * (1
    - phi));
    } */

    G::apply(models.rho_S_model, aux_data, state, out);
    G::apply(models.grav_model, state, out, tmp, cd);
    G::apply(models.mu_L_model, aux_data, out);
    G::apply(models.transport_poro_model, aux_data, tmp, state, prev_state);
    G::apply(models.perm_model, aux_data, state, out, cd, tmp);
    G::apply(models.th_osmosis_model, aux_data, out, cd);
    G::apply(models.darcy_model, aux_data, out, tmp, cd);
    G::apply(models.heat_storage_and_flux_model, aux_data, out, state, tmp, cd);
    G::apply(models.vapor_diffusion_model, aux_data, out, state, tmp, cd);
    G::apply(models.f_therm_exp_model, aux_data, tmp, state, out);
    G::apply(models.storage_model, aux_data, tmp, state, out, prev_state, cd);
    G::apply(models.eq_p_model, aux_data, state, tmp, out, cd);
    G::apply(models.eq_T_model, cd);
}

template struct ConstitutiveSetting<2>;
template struct ConstitutiveSetting<3>;

}  // namespace ConstitutiveStress_StrainTemperature
}  // namespace ProcessLib::ThermoRichardsMechanics
