/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include "MathLib/KelvinVector.h"

namespace MaterialLib
{
namespace Solids
{
namespace Phasefield
{

template <int DisplacementDim>
std::tuple<MathLib::KelvinVector::KelvinVectorType<
               DisplacementDim> /* eps_tensile */,
           MathLib::KelvinVector::KelvinVectorType<
               DisplacementDim> /* sigma_real */,
           MathLib::KelvinVector::KelvinVectorType<
               DisplacementDim> /* sigma_tensile */,
           MathLib::KelvinVector::KelvinVectorType<
               DisplacementDim> /* sigma_compressive */,
           MathLib::KelvinVector::KelvinMatrixType<DisplacementDim> /* D */,
           double /* strain_energy_tensile */, double /* elastic_energy */,
           MathLib::KelvinVector::KelvinMatrixType<
               DisplacementDim> /* C_tensile */,
           MathLib::KelvinVector::KelvinMatrixType<
               DisplacementDim> /* C_compressive
                                 */
           >
calculateOrthoVolDevDegradedStress(
    double const degradation,
    MathLib::KelvinVector::KelvinVectorType<DisplacementDim> const& eps,
    MathLib::KelvinVector::KelvinMatrixType<DisplacementDim> const& C_ortho)
{
    static constexpr int KelvinVectorSize =
        MathLib::KelvinVector::kelvin_vector_dimensions(DisplacementDim);
    using KelvinVector =
        MathLib::KelvinVector::KelvinVectorType<DisplacementDim>;
    using KelvinMatrix =
        MathLib::KelvinVector::KelvinMatrixType<DisplacementDim>;
    using Invariants = MathLib::KelvinVector::Invariants<KelvinVectorSize>;
    // calculation of deviatoric parts
    auto const& P_dev = Invariants::deviatoric_projection;

    KelvinMatrix const C = C_ortho;
    // calculation of square root of elasticity tensor
    Eigen::SelfAdjointEigenSolver<KelvinMatrix> es(C);
    KelvinMatrix const sqrt_C = es.operatorSqrt();
    Eigen::SelfAdjointEigenSolver<KelvinMatrix> es_inverse(C.inverse());
    KelvinMatrix const sqrt_C_inverse = es_inverse.operatorSqrt();

    // strain in a transformed space
    KelvinVector const epst = sqrt_C * eps;
    double const epst_curr_trace = Invariants::trace(epst);

    // projection tensors in transformed space
    KelvinMatrix teps_p = KelvinMatrix::Zero();
    KelvinMatrix teps_n = KelvinMatrix::Zero();
    if (epst_curr_trace >= 0) /* QQQ */
    {
        teps_p.template topLeftCorner<3, 3>().setConstant(1. / 3);
    }
    else
    {
        teps_n.template topLeftCorner<3, 3>().setConstant(1. / 3);
    }

    teps_p.noalias() += P_dev * KelvinMatrix::Identity();

    // strain tensile / compressive
    KelvinVector const eps_tensile = sqrt_C_inverse * (teps_p * sqrt_C) * eps;
    KelvinVector const eps_compressive =
        sqrt_C_inverse * (teps_n * sqrt_C) * eps;

    // projection tensors in original space
    KelvinMatrix const der_eps_p = sqrt_C_inverse * (teps_p * sqrt_C);
    KelvinMatrix const der_eps_n = sqrt_C_inverse * (teps_n * sqrt_C);
    // C tensile / compressive
    KelvinMatrix const C_tensile = der_eps_p.transpose() * C * der_eps_p;
    KelvinMatrix const C_compressive = der_eps_n.transpose() * C * der_eps_n;

    // stress tensile / compressive
    KelvinVector const sigma_tensile = C_tensile * eps;
    KelvinVector const sigma_compressive = C_compressive * eps;

    // decomposition of strain energy density
    double const strain_energy_tensile =
        0.5 * sigma_tensile.adjoint() * eps_tensile;
    double const strain_energy_compressive =
        0.5 * sigma_compressive.adjoint() * eps_compressive;

    double const elastic_energy =
        degradation * strain_energy_tensile + strain_energy_compressive;

    KelvinVector const sigma_real =
        degradation * sigma_tensile + sigma_compressive;
    KelvinMatrix const D = degradation * C_tensile + C_compressive;

    return std::make_tuple(eps_tensile, sigma_real, sigma_tensile,
                           sigma_compressive, D, strain_energy_tensile,
                           elastic_energy, C_tensile, C_compressive);
}

template <int DisplacementDim>
std::tuple<MathLib::KelvinVector::KelvinVectorType<
               DisplacementDim> /* eps_tensile */,
           MathLib::KelvinVector::KelvinVectorType<
               DisplacementDim> /* sigma_real */,
           MathLib::KelvinVector::KelvinVectorType<
               DisplacementDim> /* sigma_tensile */,
           MathLib::KelvinVector::KelvinMatrixType<DisplacementDim> /* D */,
           double /* strain_energy_tensile */, double /* elastic_energy */,
           MathLib::KelvinVector::KelvinMatrixType<
               DisplacementDim> /* C_tensile */,
           MathLib::KelvinVector::KelvinMatrixType<
               DisplacementDim> /* C_compressive
                                 */
           >
calculateOrthoMasonryDegradedStress(
    double const degradation,
    MathLib::KelvinVector::KelvinVectorType<DisplacementDim> const& eps,
    MathLib::KelvinVector::KelvinMatrixType<DisplacementDim> const& C_ortho)
{
    using KelvinVector =
        MathLib::KelvinVector::KelvinVectorType<DisplacementDim>;
    using KelvinMatrix =
        MathLib::KelvinVector::KelvinMatrixType<DisplacementDim>;

    KelvinMatrix const C = C_ortho;

    // calculation of square root of elasticity tensor
    Eigen::SelfAdjointEigenSolver<KelvinMatrix> es(C);
    KelvinMatrix const sqrt_C = es.operatorSqrt();
    Eigen::SelfAdjointEigenSolver<KelvinMatrix> es_inverse(C.inverse());
    KelvinMatrix const sqrt_C_inverse = es_inverse.operatorSqrt();

    // strain in a transformed space
    KelvinVector const epst = sqrt_C * eps;

    // strain tensor in 3D
    Eigen::Matrix3d const eps_3D =
        MathLib::KelvinVector::kelvinVectorToTensor(epst);

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> const es_eps_3D(eps_3D);
    Eigen::Vector3d const eigen_values_eps_3D = es_eps_3D.eigenvalues().real();
    Eigen::Matrix3d const eigen_vectors_eps_3D = es_eps_3D.eigenvectors();
    Eigen::Matrix3d const E1_eigenp =
        eigen_vectors_eps_3D.col(0) * eigen_vectors_eps_3D.col(0).transpose();

    Eigen::Matrix3d const E3_eigenp =
        eigen_vectors_eps_3D.col(2) * eigen_vectors_eps_3D.col(2).transpose();

    KelvinVector const E1_vec =
        MathLib::KelvinVector::tensorToKelvin<DisplacementDim>(E1_eigenp);

    KelvinVector const E3_vec =
        MathLib::KelvinVector::tensorToKelvin<DisplacementDim>(E3_eigenp);

    KelvinMatrix const E1oE1 = E1_vec * E1_vec.transpose();
    KelvinMatrix const E3oE3 = E3_vec * E3_vec.transpose();

    KelvinMatrix I_p = KelvinMatrix::Zero();
    KelvinMatrix I_n = KelvinMatrix::Zero();

    KelvinMatrix const I_S = KelvinMatrix::Identity();
    if (DisplacementDim == 2)
    {
        if (std::abs(eigen_values_eps_3D(2) - eigen_values_eps_3D(0)) >
            std::numeric_limits<double>::epsilon())
        {
            I_p = (macaulayTensile(eigen_values_eps_3D(0)) -
                   macaulayTensile(eigen_values_eps_3D(2))) /
                      (eigen_values_eps_3D(0) - eigen_values_eps_3D(2)) *
                      (I_S - (E1oE1 + E3oE3)) +
                  (heaviside(eigen_values_eps_3D(0)) * E1oE1 +
                   heaviside(eigen_values_eps_3D(2)) * E3oE3);
            I_n = (macaulayCompressive(eigen_values_eps_3D(0)) -
                   macaulayCompressive(eigen_values_eps_3D(2))) /
                      (eigen_values_eps_3D(0) - eigen_values_eps_3D(2)) *
                      (I_S - (E1oE1 + E3oE3)) +
                  (heaviside(-eigen_values_eps_3D(0)) * E1oE1 +
                   heaviside(-eigen_values_eps_3D(2)) * E3oE3);
        }
        else
        {
            I_p = heaviside(eigen_values_eps_3D(0)) * I_S;
            I_n = heaviside(-eigen_values_eps_3D(0)) * I_S;
        }
    }

    // strain tensile / compressive
    KelvinVector const eps_tensile = sqrt_C_inverse * (I_p * sqrt_C) * eps;
    KelvinVector const eps_compressive = sqrt_C_inverse * (I_n * sqrt_C) * eps;

    // C tensile / compressive
    KelvinMatrix const C_tensile = sqrt_C * I_p * sqrt_C;
    KelvinMatrix const C_compressive = sqrt_C * I_n * sqrt_C;

    // stress tensile / compressive
    KelvinVector const sigma_tensile = C * eps_tensile;
    KelvinVector const sigma_compressive = C * eps_compressive;

    // decomposition of strain energy density
    double const strain_energy_tensile =
        0.5 * sigma_tensile.adjoint() * eps_tensile;
    double const strain_energy_compressive =
        0.5 * sigma_compressive.adjoint() * eps_compressive;

    double const elastic_energy =
        degradation * strain_energy_tensile + strain_energy_compressive;

    KelvinVector const sigma_real =
        degradation * sigma_tensile + sigma_compressive;

    KelvinMatrix const D = degradation * C_tensile + C_compressive;

    return std::make_tuple(eps_tensile, sigma_real, sigma_tensile, D,
                           strain_energy_tensile, elastic_energy, C_tensile,
                           C_compressive);
}

}  // namespace Phasefield
}  // namespace Solids
}  // namespace MaterialLib
