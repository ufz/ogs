// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "ConstitutiveSetting.h"

#include "BaseLib/Error.h"
#include "ProcessLib/Graph/Apply.h"
#include "ProcessLib/Graph/CheckEvalOrderRT.h"

namespace ProcessLib::LargeDeformation
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

    using Inputs = mp_list<SpaceTimeData, MediaData, Temperature,
                           DeformationGradientData<D>,
                           PrevState<DeformationGradientData<D>>>;

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
    DeformationGradientData<DisplacementDim> const& deformation_gradient_data,
    GradientVectorType const& deformation_gradient_prev,
    StatefulData<DisplacementDim>& state,
    StatefulDataPrev<DisplacementDim> const& prev_state,
    MaterialStateData<DisplacementDim>& mat_state,
    ConstitutiveTempData<DisplacementDim>& tmp,
    ConstitutiveData<DisplacementDim>& cd) const
{
    auto& deformation_gradient_data_prev =
        std::get<PrevState<DeformationGradientData<DisplacementDim>>>(tmp);
    deformation_gradient_data_prev->deformation_gradient =
        deformation_gradient_prev;

    auto const aux_data =
        std::tuple{SpaceTimeData{x_position, t, dt}, MediaData{medium},
                   Temperature{T_ref}, deformation_gradient_data};
    auto const mat_state_tuple = std::tie(mat_state);

    ProcessLib::Graph::evalAllInOrder(models, aux_data, cd, mat_state_tuple,
                                      prev_state, state, tmp);
}

template struct ConstitutiveSetting<2>;
template struct ConstitutiveSetting<3>;
}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::LargeDeformation
