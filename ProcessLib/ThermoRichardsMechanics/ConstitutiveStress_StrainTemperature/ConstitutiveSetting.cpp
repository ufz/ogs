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

#include <Eigen/LU>

#include "MaterialLib/MPL/PropertyType.h"
#include "ProcessLib/Graph/Apply.h"

namespace ProcessLib::ThermoRichardsMechanics
{
namespace ConstitutiveStress_StrainTemperature
{
template <int DisplacementDim>
void ConstitutiveSetting<DisplacementDim>::init(
    ConstitutiveModels<DisplacementDim>& models, double const t,
    double const dt, ParameterLib::SpatialPosition const& x_position,
    MediaData const& media_data, TemperatureData<DisplacementDim> const& T_data,
    StatefulData<DisplacementDim>& state,
    StatefulDataPrev<DisplacementDim>& prev_state) const
{
    // Set eps_m_prev from potentially non-zero eps and sigma_sw from
    // restart.
    SpaceTimeData const x_t{x_position, t, dt};
    ElasticTangentStiffnessData<DisplacementDim> C_el_data;
    std::get<ElasticTangentStiffnessModel<DisplacementDim>>(models).eval(
        x_t, T_data, C_el_data);

    auto const& eps = std::get<StrainData<DisplacementDim>>(state).eps;
    auto const& sigma_sw =
        std::get<SwellingDataStateful<DisplacementDim>>(state).sigma_sw;
    std::get<PrevState<MechanicalStrainData<DisplacementDim>>>(prev_state)
        ->eps_m.noalias() =
        media_data.solid.hasProperty(
            MaterialPropertyLib::PropertyType::swelling_stress_rate)
            ? eps + C_el_data.C_el.inverse() * sigma_sw
            : eps;
}

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

    G::eval(std::get<ElasticTangentStiffnessModel<DisplacementDim>>(models),
            aux_data, tmp);
    G::eval(std::get<BiotModel>(models), aux_data, tmp);
    // TODO check if model eval order is still correct under all circumstances!
    G::eval(std::get<SolidCompressibilityModel<
                DisplacementDim, SolidConstitutiveRelation<DisplacementDim>>>(
                models),
            aux_data, tmp, cd);
    G::eval(std::get<SaturationModel<DisplacementDim>>(models), aux_data, state,
            tmp);

    G::eval(std::get<BishopsModel>(models), aux_data, state, tmp);
    // TODO why not ordinary state tracking?
    G::eval(std::get<BishopsPrevModel>(models), aux_data, prev_state, tmp);
    G::eval(std::get<PorosityModel<DisplacementDim>>(models), aux_data, tmp,
            state, prev_state);

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

    G::eval(std::get<SwellingModel<DisplacementDim>>(models), aux_data, state,
            prev_state, tmp);
    G::eval(std::get<SolidThermalExpansionModel<DisplacementDim>>(models),
            aux_data, tmp);
    G::eval(std::get<SolidMechanicsModel<DisplacementDim>>(models), aux_data,
            tmp, state, prev_state, mat_state_tuple, cd);

    G::eval(std::get<LiquidDensityModel<DisplacementDim>>(models), aux_data,
            out);

    /* {
        double const p_FR = -bishops_data.chi_S_L * p_cap_data.p_cap;
        // p_SR
        // TODO used by no MPL model
        variables.solid_grain_pressure =
            p_FR -
            Invariants::trace(std::get<s_mech_data>(state).sigma_eff) / (3 * (1
    - phi));
    } */

    G::eval(std::get<SolidDensityModel<DisplacementDim>>(models), aux_data,
            state, out);
    G::eval(std::get<GravityModel<DisplacementDim>>(models), state, out, tmp,
            cd);
    G::eval(std::get<LiquidViscosityModel<DisplacementDim>>(models), aux_data,
            out);
    G::eval(std::get<TransportPorosityModel<DisplacementDim>>(models), aux_data,
            tmp, state, prev_state);
    G::eval(std::get<PermeabilityModel<DisplacementDim>>(models), aux_data,
            state, out, cd, tmp);
    G::eval(std::get<ThermoOsmosisModel<DisplacementDim>>(models), aux_data,
            out, cd);
    G::eval(std::get<DarcyLawModel<DisplacementDim>>(models), aux_data, out,
            tmp, cd);
    G::eval(std::get<TRMHeatStorageAndFluxModel<DisplacementDim>>(models),
            aux_data, out, state, tmp, cd);
    G::eval(std::get<TRMVaporDiffusionModel<DisplacementDim>>(models), aux_data,
            out, state, tmp, cd);
    G::eval(std::get<FluidThermalExpansionModel<DisplacementDim>>(models),
            aux_data, tmp, state, out);
    G::eval(std::get<TRMStorageModel<DisplacementDim>>(models), aux_data, tmp,
            state, out, prev_state, cd);
    G::eval(std::get<EqPModel<DisplacementDim>>(models), aux_data, state, tmp,
            out, cd);
    G::eval(std::get<EqTModel<DisplacementDim>>(models), cd);
}

template struct ConstitutiveSetting<2>;
template struct ConstitutiveSetting<3>;

}  // namespace ConstitutiveStress_StrainTemperature
}  // namespace ProcessLib::ThermoRichardsMechanics
