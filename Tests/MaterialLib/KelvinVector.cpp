/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include <gtest/gtest.h>

#include "MaterialLib/SolidModels/KelvinVector.h"

#include "Tests/MathLib/AutoCheckTools.h"

using namespace MaterialLib::SolidModels;
namespace ac = autocheck;

template <int Size>
using KelvinVector = Eigen::Matrix<double, Size, 1, Eigen::ColMajor, Size, 1>;

template <int Size>
Eigen::Matrix<double, 3, 3> kelvinToTensor(KelvinVector<Size> const& v);

template <>
Eigen::Matrix<double, 3, 3> kelvinToTensor(KelvinVector<4> const& v)
{
    Eigen::Matrix<double, 3, 3> m;
    m << v[0], v[3] / std::sqrt(2.), 0, v[3] / std::sqrt(2.), v[1], 0, 0, 0,
        v[2];
    return m;
}

template <>
Eigen::Matrix<double, 3, 3> kelvinToTensor(KelvinVector<6> const& v)
{
    Eigen::Matrix<double, 3, 3> m;
    m << v[0], v[3] / std::sqrt(2.), v[5] / std::sqrt(2.), v[3] / std::sqrt(2.),
        v[1], v[4] / std::sqrt(2.), v[5] / std::sqrt(2.), v[4] / std::sqrt(2.),
        v[2];
    return m;
}

template <int Size>
KelvinVector<Size> tensorToKelvin(Eigen::Matrix<double, 3, 3> const& m);

template <>
KelvinVector<4> tensorToKelvin(Eigen::Matrix<double, 3, 3> const& m)
{
    EXPECT_NEAR(m(0, 1), m(1, 0), std::numeric_limits<double>::epsilon());
    EXPECT_EQ(m(0, 2), 0);
    EXPECT_EQ(m(1, 2), 0);
    EXPECT_EQ(m(2, 0), 0);
    EXPECT_EQ(m(2, 1), 0);

    KelvinVector<4> v;
    v << m(0, 0), m(1, 1), m(2, 2), m(0, 1) * std::sqrt(2.);
    return v;
}

template <>
KelvinVector<6> tensorToKelvin(Eigen::Matrix<double, 3, 3> const& m)
{
    EXPECT_NEAR(m(0, 1), m(1, 0), std::numeric_limits<double>::epsilon());
    EXPECT_NEAR(m(1, 2), m(2, 1), std::numeric_limits<double>::epsilon());
    EXPECT_NEAR(m(0, 2), m(2, 0), std::numeric_limits<double>::epsilon());

    KelvinVector<6> v;
    v << m(0, 0), m(1, 1), m(2, 2), m(0, 1) * std::sqrt(2.),
        m(1, 2) * std::sqrt(2.), m(0, 2) * std::sqrt(2.);
    return v;
}

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
    auto f = [](KelvinVector<4> const& v) {
        return (v - tensorToKelvin<4>(kelvinToTensor(v))).norm() <=
               2 * std::numeric_limits<double>::epsilon() * v.norm();
    };

    ac::check<KelvinVector<4>>(
        f, 1000, ac::make_arbitrary(kelvinVectorGenerator), gtest_reporter);
}

TEST_F(MaterialLibSolidsKelvinVector6, SelfTestMappingKelvinToTensor)
{
    auto f = [](KelvinVector<6> const& v) {
        return (v - tensorToKelvin<6>(kelvinToTensor(v))).norm() <=
               1.5 * std::numeric_limits<double>::epsilon() * v.norm();
    };

    ac::check<KelvinVector<6>>(
        f, 1000, ac::make_arbitrary(kelvinVectorGenerator), gtest_reporter);
}

//
// Determinants of a Kelvin vector and corresponding Eigen::Matrix are equal.
//

TEST_F(MaterialLibSolidsKelvinVector4, Determinant)
{
    auto f = [](KelvinVector<4> const& v) {
        return std::abs(Invariants<4>::determinant(v) -
                        kelvinToTensor(v).determinant()) <=
               std::numeric_limits<double>::epsilon() *
                   std::pow(v.norm(), 3.07);
    };

    ac::check<KelvinVector<4>>(
        f, 1000, ac::make_arbitrary(kelvinVectorGenerator), gtest_reporter);
}

TEST_F(MaterialLibSolidsKelvinVector6, Determinant)
{
    auto f = [](KelvinVector<6> const& v) {
        return std::abs(Invariants<6>::determinant(v) -
                        kelvinToTensor(v).determinant()) <=
               std::numeric_limits<double>::epsilon() *
                   std::pow(v.norm(), 3.07);
    };

    ac::check<KelvinVector<6>>(
        f, 1000, ac::make_arbitrary(kelvinVectorGenerator), gtest_reporter);
}

//
// Inverse of a Kelvin vector and corresponding Eigen::Matrix are equal.
//

TEST_F(MaterialLibSolidsKelvinVector4, Inverse)
{
    auto f = [](KelvinVector<4> const& v) {
        auto const error =
            (inverse(v) - tensorToKelvin<4>(kelvinToTensor(v).inverse()))
                .norm();
        // The error is only weekly depending on the input vector norm.
        return error < 1e-7 && 1e-9 * std::pow(v.norm(), 1.2);
    };

    auto eps = std::numeric_limits<double>::epsilon();
    ac::check<KelvinVector<4>>(
        f, 1000,
        ac::make_arbitrary(kelvinVectorGenerator)
            .discard_if([&eps](KelvinVector<4> const& v) {
                // only invertable matrices
                return (std::abs(kelvinToTensor(v).determinant()) == 0);
            }),
        gtest_reporter);
}

TEST_F(MaterialLibSolidsKelvinVector6, Inverse)
{
    auto f = [](KelvinVector<6> const& v) {
        auto const error =
            (inverse(v) - tensorToKelvin<6>(kelvinToTensor(v).inverse()))
                .norm();
        // The error is only weekly depending on the input vector norm.
        return error < 1e-7 && 1e-9 * std::pow(v.norm(), 1.2);
    };

    auto eps = std::numeric_limits<double>::epsilon();
    ac::check<KelvinVector<6>>(
        f, 1000,
        ac::make_arbitrary(kelvinVectorGenerator)
            .discard_if([&eps](KelvinVector<6> const& v) {
                // only invertable matrices
                return (std::abs(kelvinToTensor(v).determinant()) == 0);
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
