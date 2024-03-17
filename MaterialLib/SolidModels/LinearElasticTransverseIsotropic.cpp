/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on November 6, 2023, 1:42 PM
 */

#include "LinearElasticTransverseIsotropic.h"

#include <Eigen/LU>

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
LinearElasticTransverseIsotropic<DisplacementDim>::integrateStress(
    MaterialPropertyLib::VariableArray const& variable_array_prev,
    MaterialPropertyLib::VariableArray const& variable_array, double const t,
    ParameterLib::SpatialPosition const& x, double const /*dt*/,
    typename MechanicsBase<DisplacementDim>::
        MaterialStateVariables const& /*material_state_variables*/) const
{
    auto const& eps_m = std::get<MPL::SymmetricTensor<DisplacementDim>>(
        variable_array.mechanical_strain);
    auto const& eps_m_prev = std::get<MPL::SymmetricTensor<DisplacementDim>>(
        variable_array_prev.mechanical_strain);
    auto const& sigma_prev = std::get<MPL::SymmetricTensor<DisplacementDim>>(
        variable_array_prev.stress);
    auto const& T = variable_array_prev.temperature;

    KelvinMatrix C = getElasticTensor(t, x, T);

    KelvinVector sigma = sigma_prev + C * (eps_m - eps_m_prev);

    return {std::make_tuple(
        sigma,
        std::make_unique<
            typename MechanicsBase<DisplacementDim>::MaterialStateVariables>(),
        C)};
}

template <int DisplacementDim>
typename LinearElasticTransverseIsotropic<DisplacementDim>::KelvinMatrix
LinearElasticTransverseIsotropic<DisplacementDim>::
    getElasticTensorLeftTopCorner(double const t,
                                  ParameterLib::SpatialPosition const& x) const
{
    using namespace MathLib::KelvinVector;
    double const E_i = E_i_p_(t, x)[0];
    double const E_a = E_a_p_(t, x)[0];
    double const nu_i = nu_ii_p_(t, x)[0];
    double const nu_ia = nu_ia_p_(t, x)[0];

    double const nu_ai = nu_ia * (E_a / E_i);

    double const nu_ia_p_nu_ai = nu_ia * nu_ai;
    double const nu_i_p2 = nu_i * nu_i;
    double const E_i_p2 = E_i * E_i;
    // 1/D
    double const one_over_D =
        (E_i_p2 * E_a) / (1.0 - nu_i_p2 - 2.0 * (1.0 + nu_i) * nu_ia_p_nu_ai);

    double const fac1 = one_over_D / (E_i * E_a);
    double const fac2 = one_over_D / E_i_p2;
    double const a_ii = (1.0 - nu_ia_p_nu_ai) * fac1;
    double const a_ai = (1.0 - nu_i_p2) * fac2;
    double const b_ii = (nu_i + nu_ia_p_nu_ai) * fac1;
    double const b_ai = (nu_ia * (1.0 + nu_i)) * fac2;

    KelvinMatrixType<DisplacementDim> C_ortho =
        KelvinMatrixType<DisplacementDim>::Zero();

    if (DisplacementDim == 2)
    {
        // The following data are set based on the condition that the direction
        // of anisotropy is parallel to the second base (base2) of the 2D local
        // coordinate system.
        C_ortho.template topLeftCorner<3, 3>() << a_ii, b_ai, b_ii, b_ai, a_ai,
            b_ai, b_ii, b_ai, a_ii;

        return C_ortho;
    }
    C_ortho.template topLeftCorner<3, 3>() << a_ii, b_ii, b_ai, b_ii, a_ii,
        b_ai, b_ai, b_ai, a_ai;

    return C_ortho;
}

template <>
typename LinearElasticTransverseIsotropic<2>::KelvinMatrix
LinearElasticTransverseIsotropic<2>::getElasticTensor(
    double const t, ParameterLib::SpatialPosition const& x,
    double const /*T*/) const
{
    using namespace MathLib::KelvinVector;

    double const G_a = G_ia_p_(t, x)[0];
    double const c_ai = G_a;

    KelvinMatrixType<2> C_ortho = getElasticTensorLeftTopCorner(t, x);

    C_ortho.template bottomRightCorner<1, 1>().diagonal() << 2 * c_ai;

    auto const Q = [this, &x]() -> KelvinMatrixType<2>
    {
        if (!local_coordinate_system_)
        {
            return KelvinMatrixType<2>::Identity();
        }

        Eigen::Matrix<double, 2, 2> R = Eigen::Matrix<double, 2, 2>::Identity();
        R.template topLeftCorner<2, 2>().noalias() =
            local_coordinate_system_->transformation<2>(x);
        return fourthOrderRotationMatrix(R);
    }();

    // Rotate the forth-order tenser in Kelvin mapping with Q*C_ortho*Q^T and
    // return the top left corner block of size 4x4 for two-dimensional case.
    return Q * C_ortho * Q.transpose();
}

template <>
typename LinearElasticTransverseIsotropic<3>::KelvinMatrix
LinearElasticTransverseIsotropic<3>::getElasticTensor(
    double const t, ParameterLib::SpatialPosition const& x,
    double const /*T*/) const
{
    using namespace MathLib::KelvinVector;

    double const E_i = E_i_p_(t, x)[0];
    double const nu_i = nu_ii_p_(t, x)[0];
    double const G_a = G_ia_p_(t, x)[0];

    double const c_ii = E_i / (2.0 * (1 + nu_i));
    double const c_ai = G_a;

    KelvinMatrixType<3> C_ortho = getElasticTensorLeftTopCorner(t, x);

    C_ortho.template bottomRightCorner<3, 3>().diagonal() << 2.0 * c_ii,
        2 * c_ai, 2 * c_ai;

    auto const Q = [this, &x]() -> KelvinMatrixType<3>
    {
        if (!local_coordinate_system_)
        {
            return KelvinMatrixType<3>::Identity();
        }
        Eigen::Matrix3d R = Eigen::Matrix3d::Identity();
        R.template topLeftCorner<3, 3>().noalias() =
            local_coordinate_system_->transformation<3>(x);
        return fourthOrderRotationMatrix(R);
    }();

    // Rotate the forth-order tenser in Kelvin mapping with Q*C_ortho*Q^T and
    // return the 6x6 matrix for in the three-dimensional case.
    return Q * C_ortho * Q.transpose();
}

template class LinearElasticTransverseIsotropic<2>;
template class LinearElasticTransverseIsotropic<3>;

}  // namespace Solids
}  // namespace MaterialLib
