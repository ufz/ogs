/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include <gtest/gtest.h>

#include "MathLib/KelvinVector.h"

#include "Tests/AutoCheckTools.h"

using namespace MathLib::KelvinVector;
namespace ac = autocheck;

struct MaterialLibSolidsKelvinVector4 : public ::testing::Test
{
    static const int size = 4;
    ac::randomEigenMatrixGenerator<double, size, 1> kelvinVectorGenerator;

    ac::gtest_reporter gtest_reporter;
};

struct MaterialLibSolidsKelvinVector6 : public ::testing::Test
{
    static const int size = 6;
    ac::randomEigenMatrixGenerator<double, size, 1> kelvinVectorGenerator;

    ac::gtest_reporter gtest_reporter;
};

//
// Self-test of Kelvin to tensor and back conversion, which should be identity
// up to numerical errors.
//

TEST_F(MaterialLibSolidsKelvinVector4, SelfTestMappingKelvinToTensor)
{
    auto f = [](KelvinVectorType<2> const& v) {
        return (v - tensorToKelvin<2>(kelvinVectorToTensor(v))).norm() <=
               2 * std::numeric_limits<double>::epsilon() * v.norm();
    };

    ac::check<KelvinVectorType<2>>(
        f, 1000, ac::make_arbitrary(kelvinVectorGenerator), gtest_reporter);
}

TEST_F(MaterialLibSolidsKelvinVector6, SelfTestMappingKelvinToTensor)
{
    auto f = [](KelvinVectorType<3> const& v) {
        return (v - tensorToKelvin<3>(kelvinVectorToTensor(v))).norm() <=
               1.5 * std::numeric_limits<double>::epsilon() * v.norm();
    };

    ac::check<KelvinVectorType<3>>(
        f, 1000, ac::make_arbitrary(kelvinVectorGenerator), gtest_reporter);
}

//
// Determinants of a Kelvin vector and corresponding Eigen::Matrix are equal.
//

TEST_F(MaterialLibSolidsKelvinVector4, Determinant)
{
    auto f = [](KelvinVectorType<2> const& v) {
        return std::abs(Invariants<4>::determinant(v) -
                        kelvinVectorToTensor(v).determinant()) <=
               std::numeric_limits<double>::epsilon() *
                   std::pow(v.norm(), 3.07);
    };

    ac::check<KelvinVectorType<2>>(
        f, 1000, ac::make_arbitrary(kelvinVectorGenerator), gtest_reporter);
}

TEST_F(MaterialLibSolidsKelvinVector6, Determinant)
{
    auto f = [](KelvinVectorType<3> const& v) {
        return std::abs(Invariants<6>::determinant(v) -
                        kelvinVectorToTensor(v).determinant()) <=
               std::numeric_limits<double>::epsilon() *
                   std::pow(v.norm(), 3.07);
    };

    ac::check<KelvinVectorType<3>>(
        f, 1000, ac::make_arbitrary(kelvinVectorGenerator), gtest_reporter);
}

//
// Inverse of a Kelvin vector and corresponding Eigen::Matrix are equal.
//

TEST_F(MaterialLibSolidsKelvinVector4, Inverse)
{
    auto f = [](KelvinVectorType<2> const& v) {
        auto const error =
            (inverse(v) - tensorToKelvin<2>(kelvinVectorToTensor(v).inverse()))
                .norm();
        // The error is only weekly depending on the input vector norm.
        return error < 1e-6 && error < 1e-8 * std::pow(v.norm(), 1.4);
    };

    ac::check<KelvinVectorType<2>>(
        f, 1000,
        ac::make_arbitrary(kelvinVectorGenerator)
            .discard_if([](KelvinVectorType<2> const& v) {
                // only invertable matrices
                return (std::abs(kelvinVectorToTensor(v).determinant()) == 0);
            }),
        gtest_reporter);
}

TEST_F(MaterialLibSolidsKelvinVector6, Inverse)
{
    auto f = [](KelvinVectorType<3> const& v) {
        auto const error =
            (inverse(v) - tensorToKelvin<3>(kelvinVectorToTensor(v).inverse()))
                .norm();
        // The error is only weekly depending on the input vector norm.
        return error < 1e-6 && error < 1e-8 * std::pow(v.norm(), 1.4);
    };

    ac::check<KelvinVectorType<3>>(
        f, 1000,
        ac::make_arbitrary(kelvinVectorGenerator)
            .discard_if([](KelvinVectorType<3> const& v) {
                // only invertable matrices
                return (std::abs(kelvinVectorToTensor(v).determinant()) == 0);
            }),
        gtest_reporter);
}

//
// Tests of deviatoric and spherical projection tensors.
//

TEST_F(MaterialLibSolidsKelvinVector4, DeviatoricSphericalProjections)
{
    auto const& P_dev = Invariants<size>::deviatoric_projection;
    auto const& P_sph = Invariants<size>::spherical_projection;

    // Test product P_dev * P_sph is zero.
    Eigen::Matrix<double, 4, 4> const P_dev_P_sph = P_dev * P_sph;
    EXPECT_LT(P_dev_P_sph.norm(), std::numeric_limits<double>::epsilon());
    EXPECT_LT(P_dev_P_sph.maxCoeff(), std::numeric_limits<double>::epsilon());

    // Test product P_sph * P_dev is zero.
    Eigen::Matrix<double, 4, 4> const P_sph_P_dev = P_sph * P_dev;
    EXPECT_LT(P_sph_P_dev.norm(), std::numeric_limits<double>::epsilon());
    EXPECT_LT(P_sph_P_dev.maxCoeff(), std::numeric_limits<double>::epsilon());

    // Test sum is identity.
    EXPECT_EQ(P_dev + P_sph, (Eigen::Matrix<double, size, size>::Identity()));
}

TEST_F(MaterialLibSolidsKelvinVector6, DeviatoricSphericalProjections)
{
    auto const& P_dev = Invariants<size>::deviatoric_projection;
    auto const& P_sph = Invariants<size>::spherical_projection;

    // Test product P_dev * P_sph is zero.
    Eigen::Matrix<double, 6, 6> const P_dev_P_sph = P_dev * P_sph;
    EXPECT_LT(P_dev_P_sph.norm(), std::numeric_limits<double>::epsilon());
    EXPECT_LT(P_dev_P_sph.maxCoeff(), std::numeric_limits<double>::epsilon());

    // Test product P_sph * P_dev is zero.
    Eigen::Matrix<double, 6, 6> const P_sph_P_dev = P_sph * P_dev;
    EXPECT_LT(P_sph_P_dev.norm(), std::numeric_limits<double>::epsilon());
    EXPECT_LT(P_sph_P_dev.maxCoeff(), std::numeric_limits<double>::epsilon());

    // Test sum is identity.
    EXPECT_EQ(P_dev + P_sph, (Eigen::Matrix<double, size, size>::Identity()));
}
