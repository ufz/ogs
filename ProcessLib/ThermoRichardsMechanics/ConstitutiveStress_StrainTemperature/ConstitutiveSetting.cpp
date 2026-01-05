// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "ConstitutiveSetting.h"

#include <Eigen/LU>

#include "MaterialLib/MPL/PropertyType.h"
#include "ProcessLib/Graph/Apply.h"
#include "ProcessLib/Graph/CheckEvalOrderRT.h"

namespace ProcessLib::ThermoRichardsMechanics
{
namespace ConstitutiveStress_StrainTemperature
{
template <int DisplacementDim>
static bool checkCorrectModelEvalOrder()
{
    INFO(
        "Checking correct model evaluation order in the constitutive setting.");

    using namespace boost::mp11;

    constexpr auto D = DisplacementDim;

    using Inputs = mp_list<SpaceTimeData, MediaData, TemperatureData<D>,
                           CapillaryPressureData<D>, StrainData<D>>;

    using InputsAndPrevState = mp_append<Inputs, StatefulDataPrev<D>>;

    bool const is_correct = ProcessLib::Graph::isEvalOrderCorrectRT<
        ConstitutiveModels<DisplacementDim>, InputsAndPrevState>();

    if (!is_correct)
    {
        OGS_FATAL("The constitutive setting has a wrong evaluation order.");
    }

    INFO("Model evaluation order is correct.");

    return is_correct;
}

template <int DisplacementDim>
void ConstitutiveSetting<DisplacementDim>::init(
    ConstitutiveModels<DisplacementDim>& models, double const t,
    double const dt, ParameterLib::SpatialPosition const& x_position,
    MediaData const& media_data, TemperatureData<DisplacementDim> const& T_data,
    StatefulData<DisplacementDim>& state,
    StatefulDataPrev<DisplacementDim>& prev_state) const
{
    [[maybe_unused]] static const bool model_order_correct =
        checkCorrectModelEvalOrder<DisplacementDim>();

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
    auto const aux_data = std::tuple{SpaceTimeData{x_position, t, dt},
                                     MediaData{medium}, T_data, p_cap_data};

    auto const mat_state_tuple = std::tie(mat_state);

    // TODO will eps lag one iteration behind? (since it's not updated after
    // solving the global equation system)
    std::get<StrainData<DisplacementDim>>(state).eps.noalias() = eps_arg;

    ProcessLib::Graph::evalAllInOrder(models, aux_data, cd, mat_state_tuple,
                                      out, prev_state, state, tmp);

    // TODO why not ordinary state tracking for BishopsPrevModel?

    {
        auto const& biot_data = std::get<BiotData>(tmp);
        auto const& poro_data = std::get<PorosityData>(state);

        if (biot_data() < poro_data.phi)
        {
            OGS_FATAL(
                "ThermoRichardsMechanics: Biot-coefficient {} is smaller than "
                "porosity {} in element {}.",
                biot_data(), poro_data.phi, *x_position.getElementID());
        }
    }
}

template struct ConstitutiveSetting<2>;
template struct ConstitutiveSetting<3>;

}  // namespace ConstitutiveStress_StrainTemperature
}  // namespace ProcessLib::ThermoRichardsMechanics
