/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ConstitutiveSetting.h"

#include "ProcessLib/Graph/Apply.h"
#include "ProcessLib/Graph/CheckEvalOrderRT.h"

namespace ProcessLib::ThermoRichardsMechanics
{
namespace ConstitutiveStressSaturation_StrainPressureTemperature
{
template <int DisplacementDim>
static bool checkCorrectModelEvalOrder()
{
    INFO(
        "Checking correct model evaluation order in the constitutive setting.");

    using namespace boost::mp11;

    constexpr auto D = DisplacementDim;

    using Inputs =
        mp_list<SpaceTimeData, MediaData, TemperatureData<D>,
                CapillaryPressureData<D>, StrainData<D>
                //, MaterialStateData<D> /*TODO material state data is a special
                // case: it's both input and output data.*/
                >;

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
    ConstitutiveModels<DisplacementDim>&, double const /*t*/,
    double const /*dt*/, ParameterLib::SpatialPosition const&, MediaData const&,
    TemperatureData<DisplacementDim> const&, StatefulData<DisplacementDim>&,
    StatefulDataPrev<DisplacementDim>&) const
{
    [[maybe_unused]] static const bool model_order_correct =
        checkCorrectModelEvalOrder<DisplacementDim>();
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
                "porosity {} in element/integration point {}/{}.",
                biot_data(), poro_data.phi, *x_position.getElementID(),
                *x_position.getIntegrationPoint());
        }
    }

    // TODO Solid thermal expansion is not needed for solid mechanics (it is
    // computed by the solid material model itself), but for fluid expansion.
    // This duplication should be avoided in the future.
}

template struct ConstitutiveSetting<2>;
template struct ConstitutiveSetting<3>;
}  // namespace ConstitutiveStressSaturation_StrainPressureTemperature
}  // namespace ProcessLib::ThermoRichardsMechanics
