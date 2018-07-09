/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "LinearElasticIsotropic.h"

namespace MaterialLib
{
namespace Solids
{
template <int DisplacementDim>
boost::optional<
    std::tuple<typename LinearElasticIsotropic<DisplacementDim>::KelvinVector,
               std::unique_ptr<typename MechanicsBase<
                   DisplacementDim>::MaterialStateVariables>,
               typename LinearElasticIsotropic<DisplacementDim>::KelvinMatrix>>
LinearElasticIsotropic<DisplacementDim>::integrateStress(
    double const t,
    ProcessLib::SpatialPosition const& x,
    double const /*dt*/,
    KelvinVector const& eps_prev,
    KelvinVector const& eps,
    KelvinVector const& sigma_prev,
    typename MechanicsBase<DisplacementDim>::MaterialStateVariables const&
        material_state_variables, double const /*T*/, double const /*p*/)
{
    KelvinMatrix C = KelvinMatrix::Zero();

    C.template topLeftCorner<3, 3>().setConstant(_mp.lambda(t, x));
    C.noalias() += 2 * _mp.mu(t, x) * KelvinMatrix::Identity();

    KelvinVector sigma = sigma_prev + C * (eps - eps_prev);

    return {std::make_tuple(
        sigma,
        std::unique_ptr<
            typename MechanicsBase<DisplacementDim>::MaterialStateVariables>{
            new MaterialStateVariables{
                static_cast<MaterialStateVariables const&>(
                    material_state_variables)}},
        C)};
}

template class LinearElasticIsotropic<2>;
template class LinearElasticIsotropic<3>;

}   // namespace Solids
}   // namespace MaterialLib
