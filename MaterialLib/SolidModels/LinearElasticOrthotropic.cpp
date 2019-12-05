/**
 * \file
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <Eigen/LU>

#include "LinearElasticOrthotropic.h"

namespace MaterialLib
{
namespace Solids
{
template <int DisplacementDim>
std::optional<
    std::tuple<typename MechanicsBase<DisplacementDim>::KelvinVector,
               std::unique_ptr<typename MechanicsBase<
                   DisplacementDim>::MaterialStateVariables>,
               typename MechanicsBase<DisplacementDim>::KelvinMatrix>>
LinearElasticOrthotropic<DisplacementDim>::integrateStress(
    double const t, ParameterLib::SpatialPosition const& x, double const /*dt*/,
    KelvinVector const& eps_prev, KelvinVector const& eps,
    KelvinVector const& sigma_prev,
    typename MechanicsBase<DisplacementDim>::
        MaterialStateVariables const& /* material_state_variables */,
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
typename LinearElasticOrthotropic<DisplacementDim>::KelvinMatrix
LinearElasticOrthotropic<DisplacementDim>::getElasticTensor(
    double const t, ParameterLib::SpatialPosition const& x,
    double const /*T*/) const
{
    using namespace MathLib::KelvinVector;

    auto const& mp = _mp.evaluate(t, x);
    auto const E = [&mp](int const i) { return mp.E(i); };
    auto const G = [&mp](int const i, int const j) { return mp.G(i, j); };
    auto const nu = [&mp](int const i, int const j) { return mp.nu(i, j); };

    KelvinMatrixType<3> S_ortho = KelvinMatrixType<3>::Zero();
    // clang-format off
    S_ortho.template topLeftCorner<3, 3>() <<
               1. / E(1), -nu(2, 1) / E(2), -nu(3, 1) / E(3),
        -nu(1, 2) / E(1),        1. / E(2), -nu(3, 2) / E(3),
        -nu(1, 3) / E(1), -nu(2, 3) / E(2),        1. / E(3);

    S_ortho.template bottomRightCorner<3, 3>().diagonal() <<
        1. / (2 * G(1, 2)),
        1. / (2 * G(2, 3)),
        1. / (2 * G(1, 3));
    // clang-format on

    KelvinMatrixType<3> const C_ortho = S_ortho.inverse();
    auto const Q = [this, &x]() -> KelvinMatrixType<3> {
        if (!_local_coordinate_system)
        {
            return MathLib::KelvinVector::KelvinMatrixType<3>::Identity();
        }
        return MathLib::KelvinVector::fourthOrderRotationMatrix(
            _local_coordinate_system->transformation<3>(x));
    }();

    // Rotate the forth-order tenser in Kelvin mapping with Q*C_ortho*Q^T and
    // return the top left corner block of size 4x4 for two-dimensional case or
    // the full 6x6 matrix is returned in the three-dimensional case.
    return (Q * C_ortho * Q.transpose())
        .template topLeftCorner<
            KelvinVectorDimensions<DisplacementDim>::value,
            KelvinVectorDimensions<DisplacementDim>::value>();
}

template class LinearElasticOrthotropic<2>;
template class LinearElasticOrthotropic<3>;

}  // namespace Solids
}  // namespace MaterialLib
