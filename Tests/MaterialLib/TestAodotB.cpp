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

using MathLib::KelvinVector::KelvinMatrixType;
using MathLib::KelvinVector::KelvinVectorType;

namespace KV = MathLib::KelvinVector;

namespace
{

inline std::pair<int, int> pair3D(int I)
{
    switch (I)
    {
        case 0:
            return {0, 0};
        case 1:
            return {1, 1};
        case 2:
            return {2, 2};
        case 3:
            return {0, 1};
        case 4:
            return {1, 2};
        case 5:
            return {0, 2};
        default:
            return {0, 0};
    }
}

inline std::pair<int, int> pair2D(int I)
{
    switch (I)
    {
        case 0:
            return {0, 0};
        case 1:
            return {1, 1};
        case 2:
            return {2, 2};
        case 3:
            return {0, 1};
        default:
            return {0, 0};
    }
}

// T_ijkl = 1/4 (A_ik B_jl + A_il B_jk + A_jk B_il + A_jl B_ik)

inline double T_ref(Eigen::Matrix3d const& A, Eigen::Matrix3d const& B, int i,
                    int j, int k, int l)
{
    return 0.25 * (A(i, k) * B(j, l) + A(i, l) * B(j, k) + A(j, k) * B(i, l) +
                   A(j, l) * B(i, k));
}

}  // namespace

TEST(AOdotBTest, ThreeD_MatchesReference)
{
    Eigen::Matrix3d A;
    A << 1.2, 0.3, -0.4, 0.3, 2.1, 0.5, -0.4, 0.5, 0.8;
    Eigen::Matrix3d B;
    B << -0.7, 0.2, 0.9, 0.2, 1.4, -0.1, 0.9, -0.1, 0.6;

    A = 0.5 * (A + A.transpose());
    B = 0.5 * (B + B.transpose());

    const Eigen::Matrix<double, 6, 1> alpha(
        (Eigen::Matrix<double, 6, 1>() << 1.0, 1.0, 1.0, std::sqrt(2.0),
         std::sqrt(2.0), std::sqrt(2.0))
            .finished());

    KV::KelvinMatrixType<3> M_ref;
    M_ref.setZero();
    for (int I = 0; I < 6; ++I)
    {
        auto [i, j] = pair3D(I);
        for (int J = 0; J < 6; ++J)
        {
            auto [k, l] = pair3D(J);
            M_ref(I, J) = alpha(I) * alpha(J) * T_ref(A, B, i, j, k, l);
        }
    }

    KV::KelvinVectorType<3> Ak, Bk;
    Ak << A(0, 0), A(1, 1), A(2, 2), std::sqrt(2.0) * A(0, 1),
        std::sqrt(2.0) * A(1, 2), std::sqrt(2.0) * A(0, 2);
    Bk << B(0, 0), B(1, 1), B(2, 2), std::sqrt(2.0) * B(0, 1),
        std::sqrt(2.0) * B(1, 2), std::sqrt(2.0) * B(0, 2);

    const auto M = MaterialLib::Solids::Phasefield::aOdotB<3>(Ak, Bk);

    const double tol = 1e-12;
    EXPECT_TRUE(((M - M_ref).array().abs() < tol).all());
}

TEST(AOdotBTest, TwoD_MatchesReference)
{
    Eigen::Matrix3d A = Eigen::Matrix3d::Zero();
    A(0, 0) = 1.0;
    A(1, 1) = 2.0;
    A(2, 2) = 0.5;
    A(0, 1) = A(1, 0) = 0.3;
    Eigen::Matrix3d B = Eigen::Matrix3d::Zero();
    B(0, 0) = -0.4;
    B(1, 1) = 1.1;
    B(2, 2) = 0.2;
    B(0, 1) = B(1, 0) = -0.6;

    const Eigen::Matrix<double, 4, 1> alpha2(
        (Eigen::Matrix<double, 4, 1>() << 1.0, 1.0, 1.0, std::sqrt(2.0))
            .finished());

    Eigen::Matrix<double, 4, 4> M_ref;
    M_ref.setZero();
    for (int I = 0; I < 4; ++I)
    {
        auto [i, j] = pair2D(I);
        for (int J = 0; J < 4; ++J)
        {
            auto [k, l] = pair2D(J);
            M_ref(I, J) = alpha2(I) * alpha2(J) * T_ref(A, B, i, j, k, l);
        }
    }

    KV::KelvinVectorType<2> Ak, Bk;
    Ak << A(0, 0), A(1, 1), A(2, 2), std::sqrt(2.0) * A(0, 1);
    Bk << B(0, 0), B(1, 1), B(2, 2), std::sqrt(2.0) * B(0, 1);

    const auto M = MaterialLib::Solids::Phasefield::aOdotB<2>(Ak, Bk);

    const double tol = 1e-12;
    ASSERT_EQ(M.rows(), 4);
    ASSERT_EQ(M.cols(), 4);
    for (int r = 0; r < 4; ++r)
        for (int c = 0; c < 4; ++c)
            EXPECT_NEAR(M(r, c), M_ref(r, c), tol);
}
