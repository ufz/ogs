// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "ConstitutiveSetting.h"

#include "ProcessLib/Graph/Apply.h"
#include "ProcessLib/Graph/CheckEvalOrderRT.h"

namespace ProcessLib::SmallDeformation
{
namespace ConstitutiveRelations
{
template <int DisplacementDim>
static bool checkCorrectModelEvalOrder()
{
    INFO(
        "Checking correct model evaluation order in the constitutive setting.");

    using namespace boost::mp11;

    constexpr auto D = DisplacementDim;

    using Inputs = mp_list<SpaceTimeData, MediaData, Temperature, StrainData<D>,
                           PrevState<StrainData<D>>>;

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
void ConstitutiveSetting<DisplacementDim>::init()
{
    [[maybe_unused]] static const bool model_order_correct =
        checkCorrectModelEvalOrder<DisplacementDim>();
}
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
    auto& eps_data = std::get<StrainData<DisplacementDim>>(out);
    eps_data.eps = eps;
    auto& eps_data_prev = std::get<PrevState<StrainData<DisplacementDim>>>(tmp);
    eps_data_prev->eps = eps_prev;

    auto const aux_data = std::tuple{SpaceTimeData{x_position, t, dt},
                                     MediaData{medium}, Temperature{T_ref}};
    auto const mat_state_tuple = std::tie(mat_state);

    ProcessLib::Graph::evalAllInOrder(models, aux_data, cd, mat_state_tuple,
                                      out, prev_state, state, tmp);
}

template struct ConstitutiveSetting<2>;
template struct ConstitutiveSetting<3>;
}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::SmallDeformation
