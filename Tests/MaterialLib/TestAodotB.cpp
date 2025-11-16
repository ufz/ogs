/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 * Created on November 1, 2016, 11:06 AM
 */

#include <MaterialLib/SolidModels/LinearElasticIsotropicPhaseField.h>
#include <gtest/gtest.h>

#include <Eigen/Dense>
#include <cmath>
#include <range/v3/all.hpp>

using MathLib::KelvinVector::KelvinMatrixType;
using MathLib::KelvinVector::KelvinVectorType;

namespace KV = MathLib::KelvinVector;

namespace
{

constexpr static std::pair<int, int> pairs[6] = {{0, 0}, {1, 1}, {2, 2},
                                                 {0, 1}, {1, 2}, {0, 2}};

// T_ijkl = 1/4 (A_ik B_jl + A_il B_jk + A_jk B_il + A_jl B_ik)
double T_ref(const Eigen::Matrix3d& A, const Eigen::Matrix3d& B, const int i,
             const int j, const int k, const int l)
{
    return 0.25 * (A(i, k) * B(j, l) + A(i, l) * B(j, k) + A(j, k) * B(i, l) +
                   A(j, l) * B(i, k));
}

}  // namespace

TEST(AOdotBTest, ThreeD_MatchesReference)
{
    const Eigen::Matrix3d A =
        (Eigen::Matrix3d() << 1.2, 0.3, -0.4, 0.3, 2.1, 0.5, -0.4, 0.5, 0.8)
            .finished();

    const Eigen::Matrix3d B =
        (Eigen::Matrix3d() << -0.7, 0.2, 0.9, 0.2, 1.4, -0.1, 0.9, -0.1, 0.6)
            .finished();

    const Eigen::Matrix3d Asym = 0.5 * (A + A.transpose());
    const Eigen::Matrix3d Bsym = 0.5 * (B + B.transpose());

    const Eigen::Matrix<double, 6, 1> alpha =
        (Eigen::Matrix<double, 6, 1>() << 1.0, 1.0, 1.0, std::sqrt(2.0),
         std::sqrt(2.0), std::sqrt(2.0))
            .finished();

    KV::KelvinMatrixType<3> M_ref = KV::KelvinMatrixType<3>::Zero();

    auto pairs_indices = ranges::views::indices(0, 6);
    auto cartesian =
        ranges::views::cartesian_product(pairs_indices, pairs_indices);

    for (const auto& [I, J] : cartesian)
    {
        const auto& [i, j] = pairs[I];
        const auto& [k, l] = pairs[J];
        M_ref(I, J) = alpha(I) * alpha(J) * T_ref(Asym, Bsym, i, j, k, l);
    }

    const KV::KelvinVectorType<3> Ak =
        (KV::KelvinVectorType<3>() << Asym(0, 0), Asym(1, 1), Asym(2, 2),
         std::sqrt(2.0) * Asym(0, 1), std::sqrt(2.0) * Asym(1, 2),
         std::sqrt(2.0) * Asym(0, 2))
            .finished();
    const KV::KelvinVectorType<3> Bk =
        (KV::KelvinVectorType<3>() << Bsym(0, 0), Bsym(1, 1), Bsym(2, 2),
         std::sqrt(2.0) * Bsym(0, 1), std::sqrt(2.0) * Bsym(1, 2),
         std::sqrt(2.0) * Bsym(0, 2))
            .finished();

    const auto M = MaterialLib::Solids::Phasefield::aOdotB<3>(Ak, Bk);

    const double tol = 1e-12;
    EXPECT_TRUE(((M - M_ref).array().abs() < tol).all());
}

TEST(AOdotBTest, TwoD_MatchesReference)
{
    const Eigen::Matrix3d A =
        (Eigen::Matrix3d() << 1.0, 0.3, 0.0, 0.3, 2.0, 0.0, 0.0, 0.0, 0.5)
            .finished();

    const Eigen::Matrix3d B =
        (Eigen::Matrix3d() << -0.4, -0.6, 0.0, -0.6, 1.1, 0.0, 0.0, 0.0, 0.2)
            .finished();

    const Eigen::Matrix<double, 4, 1> alpha =
        (Eigen::Matrix<double, 4, 1>() << 1.0, 1.0, 1.0, std::sqrt(2.0))
            .finished();

    Eigen::Matrix<double, 4, 4> M_ref = Eigen::Matrix<double, 4, 4>::Zero();

    auto pairs_indices_2d = ranges::views::indices(0, 4);
    auto cartesian_2d =
        ranges::views::cartesian_product(pairs_indices_2d, pairs_indices_2d);

    for (const auto& [I, J] : cartesian_2d)
    {
        const auto& [i, j] = pairs[I];
        const auto& [k, l] = pairs[J];
        M_ref(I, J) = alpha(I) * alpha(J) * T_ref(A, B, i, j, k, l);
    }

    const KV::KelvinVectorType<2> Ak =
        (KV::KelvinVectorType<2>() << A(0, 0), A(1, 1), A(2, 2),
         std::sqrt(2.0) * A(0, 1))
            .finished();
    const KV::KelvinVectorType<2> Bk =
        (KV::KelvinVectorType<2>() << B(0, 0), B(1, 1), B(2, 2),
         std::sqrt(2.0) * B(0, 1))
            .finished();

    const auto M = MaterialLib::Solids::Phasefield::aOdotB<2>(Ak, Bk);

    const double tol = 1e-12;
    EXPECT_TRUE(((M - M_ref).array().abs() < tol).all());
}
