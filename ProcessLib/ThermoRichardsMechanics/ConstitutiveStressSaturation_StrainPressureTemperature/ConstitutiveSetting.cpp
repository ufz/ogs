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

namespace ProcessLib::ThermoRichardsMechanics
{
namespace ConstitutiveStressSaturation_StrainPressureTemperature
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
    namespace G = ProcessLib::Graph;

    auto const aux_data = std::tuple{SpaceTimeData{x_position, t, dt},
                                     MediaData{medium}, T_data, p_cap_data};

    auto const mat_state_tuple = std::tie(mat_state);

    // TODO will eps lag one iteration behind? (since it's not updated after
    // solving the global equation system)
    std::get<StrainData<DisplacementDim>>(state).eps.noalias() = eps_arg;

    G::eval(models.biot_model, aux_data, tmp);

    G::eval(models.s_mech_model, aux_data, tmp, state, prev_state,
            mat_state_tuple, cd);

    G::eval(models.solid_compressibility_model, aux_data, tmp, cd);

    G::eval(models.bishops_model, aux_data, state, tmp);
    // TODO why not ordinary state tracking?
    G::eval(models.bishops_prev_model, aux_data, prev_state, tmp);
    G::eval(models.poro_model, aux_data, tmp, state, prev_state);

    // TODO move to local assembler?
    {
        auto const& biot_data = std::get<BiotData>(tmp);
        auto const& poro_data = std::get<PorosityData>(state);

        if (biot_data() < poro_data.phi)
        {
            OGS_FATAL(
                "ThermoRichardsMechanics: Biot-coefficient {} is smaller than "
                "porosity {} in element/integration point {}/{}.",
                biot_data(), poro_data.phi, *x_position.getElementID(),
                *x_position.getIntegrationPoint());
        }
    }

    G::eval(models.rho_L_model, aux_data, out);
    G::eval(models.rho_S_model, aux_data, state, out);
    G::eval(models.grav_model, state, out, tmp, cd);
    G::eval(models.mu_L_model, aux_data, out);
    G::eval(models.transport_poro_model, aux_data, tmp, state, prev_state);
    G::eval(models.perm_model, aux_data, state, out, cd, tmp);
    G::eval(models.th_osmosis_model, aux_data, out, cd);
    G::eval(models.darcy_model, aux_data, out, tmp, cd);
    G::eval(models.heat_storage_and_flux_model, aux_data, out, state, tmp, cd);
    G::eval(models.vapor_diffusion_model, aux_data, out, state, tmp, cd);

    // TODO Not needed for solid mechanics (solid thermal expansion is computed
    // by the solid material model), but for fluid expansion. This duplication
    // should be avoided in the future.
    G::eval(models.s_therm_exp_model, aux_data, tmp);

    G::eval(models.f_therm_exp_model, aux_data, tmp, state, out);
    G::eval(models.storage_model, aux_data, tmp, state, out, prev_state, cd);
    G::eval(models.eq_p_model, aux_data, state, tmp, out, cd);
    G::eval(models.eq_T_model, cd);
}

template struct ConstitutiveSetting<2>;
template struct ConstitutiveSetting<3>;
}  // namespace ConstitutiveStressSaturation_StrainPressureTemperature
}  // namespace ProcessLib::ThermoRichardsMechanics
