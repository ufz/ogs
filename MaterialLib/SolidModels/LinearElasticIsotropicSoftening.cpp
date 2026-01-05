// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "LinearElasticIsotropicSoftening.h"

#include "MaterialLib/MPL/Utils/GetSymmetricTensor.h"

namespace MPL = MaterialPropertyLib;

namespace MaterialLib
{
namespace Solids
{
template <int DisplacementDim>
std::optional<std::tuple<typename MechanicsBase<DisplacementDim>::KelvinVector,
                         std::unique_ptr<typename MechanicsBase<
                             DisplacementDim>::MaterialStateVariables>,
                         typename MechanicsBase<DisplacementDim>::KelvinMatrix>>
LinearElasticIsotropicSoftening<DisplacementDim>::integrateStress(
    MaterialPropertyLib::VariableArray const& variable_array_prev,
    MaterialPropertyLib::VariableArray const& variable_array, double const t,
    ParameterLib::SpatialPosition const& x, double const dt,
    typename MechanicsBase<DisplacementDim>::
        MaterialStateVariables const& /*material_state_variables*/) const
{
    auto const& eps_m = std::get<MPL::SymmetricTensor<DisplacementDim>>(
        variable_array.mechanical_strain);
    auto const& eps_m_prev = std::get<MPL::SymmetricTensor<DisplacementDim>>(
        variable_array_prev.mechanical_strain);
    auto const& sigma_prev = std::get<MPL::SymmetricTensor<DisplacementDim>>(
        variable_array_prev.stress);
    auto const T = variable_array_prev.temperature;

    double const strength = _strength(t, x)[0];
    double const strength_prev = _strength(t - dt, x)[0];

    KelvinMatrix C = strength * getElasticTensor(t, x, T);
    KelvinVector sigma = sigma_prev + C * (eps_m - eps_m_prev);

    using Invariants = MathLib::KelvinVector::Invariants<KelvinVectorSize>;
    double const sig_tr = Invariants::trace(sigma);
    if (strength != strength_prev && strength_prev != 0.0 && sig_tr != 0.0)
    {
        sigma *=
            strength * Invariants::trace(sigma_prev) / (strength_prev * sig_tr);
    }

    return {std::make_tuple(
        sigma,
        std::make_unique<
            typename MechanicsBase<DisplacementDim>::MaterialStateVariables>(),
        C)};
}

template <int DisplacementDim>
typename LinearElasticIsotropicSoftening<DisplacementDim>::KelvinMatrix
LinearElasticIsotropicSoftening<DisplacementDim>::getElasticTensor(
    double const t, ParameterLib::SpatialPosition const& x,
    double const /*T*/) const
{
    using KelvinMatrix =
        MathLib::KelvinVector::KelvinMatrixType<DisplacementDim>;

    KelvinMatrix tangentStiffness = KelvinMatrix::Zero();
    tangentStiffness.template topLeftCorner<3, 3>().setConstant(
        _mp.lambda(t, x));
    tangentStiffness.noalias() += 2 * _mp.mu(t, x) * KelvinMatrix::Identity();
    return tangentStiffness;
}

template class LinearElasticIsotropicSoftening<2>;
template class LinearElasticIsotropicSoftening<3>;

}  // namespace Solids
}  // namespace MaterialLib
