/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "LinearElasticIsotropic.h"

namespace MaterialLib
{
namespace Solids
{
template <int DisplacementDim>
std::optional<std::tuple<typename MechanicsBase<DisplacementDim>::KelvinVector,
                         std::unique_ptr<typename MechanicsBase<
                             DisplacementDim>::MaterialStateVariables>,
                         typename MechanicsBase<DisplacementDim>::KelvinMatrix>>
LinearElasticIsotropic<DisplacementDim>::integrateStress(
    MaterialPropertyLib::VariableArray const& variable_array_prev,
    MaterialPropertyLib::VariableArray const& variable_array, double const t,
    ParameterLib::SpatialPosition const& x, double const /*dt*/,
    KelvinVector const& eps_prev, KelvinVector const& eps,
    KelvinVector const& sigma_prev,
    typename MechanicsBase<DisplacementDim>::
        MaterialStateVariables const& /*material_state_variables*/,
    double const T) const
{
    KelvinMatrix C = getElasticTensor(t, x, T);

    KelvinVector sigma = sigma_prev + C * (eps - eps_prev);

    return {std::make_tuple(
        sigma,
        std::make_unique<
            typename MechanicsBase<DisplacementDim>::MaterialStateVariables>(),
        C)};
}

template <int DisplacementDim>
typename LinearElasticIsotropic<DisplacementDim>::KelvinMatrix
LinearElasticIsotropic<DisplacementDim>::getElasticTensor(
    double const t, ParameterLib::SpatialPosition const& x,
    double const /*T*/) const
{
    return elasticTangentStiffness<DisplacementDim>(_mp.lambda(t, x),
                                                    _mp.mu(t, x));
}

template class LinearElasticIsotropic<2>;
template class LinearElasticIsotropic<3>;

}   // namespace Solids
}   // namespace MaterialLib
