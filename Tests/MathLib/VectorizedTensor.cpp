/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "MathLib/VectorizedTensor.h"

#include <gtest/gtest.h>

#include <Eigen/LU>
#include <range/v3/algorithm/for_each.hpp>
#include <range/v3/view/cartesian_product.hpp>
#include <range/v3/view/iota.hpp>

#include "Tests/AutoCheckTools.h"

using namespace MathLib;
namespace ac = autocheck;

TEST(VectorizedTensorTest, DimensionsAndSizes)
{
    namespace VT = MathLib::VectorizedTensor;

    // Type implicitly uses size(), which is being tested.
    static_assert(VT::dimension<VT::Type<1>>() == 1);
    static_assert(VT::dimension<VT::Type<2>>() == 2);
    static_assert(VT::dimension<VT::Type<3>>() == 3);
}

//
// Self-tests of a tensor in matrix form to a vector form and back conversion.
//
template <int DisplacementDim>
auto vectorizedTensorToMatrixMapping()
{
    static constexpr int size = VectorizedTensor::size(DisplacementDim);
    ac::randomEigenMatrixGenerator<double, size, 1> generator;

    auto to_tensor_and_back =
        [](VectorizedTensor::Type<DisplacementDim> const& v)
    {
        return (v == VectorizedTensor::toVector<DisplacementDim>(
                         VectorizedTensor::toTensor<DisplacementDim>(v)));
    };

    return ac::check<VectorizedTensor::Type<DisplacementDim>>(
        to_tensor_and_back, 1000, ac::make_arbitrary(generator),
        ac::gtest_reporter());
}

template <int DisplacementDim, typename Generator>
auto matrixToVectorizedTensorMapping(Generator const& generator)
{
    auto to_vector_and_back = [](Eigen::Matrix3d const& t)
    {
        return (t == VectorizedTensor::toTensor<DisplacementDim>(
                         VectorizedTensor::toVector<DisplacementDim>(t)));
    };

    return ac::check<Eigen::Matrix3d>(to_vector_and_back, 1000,
                                      ac::make_arbitrary(generator),
                                      ac::gtest_reporter());
}

struct random1dConvertibleEigenMatrixGenerator
{
    ac::generator<double> source;
    using result_type = Eigen::Matrix<double, 3, 3>;

    result_type operator()(std::size_t size = 0)
    {
        result_type rv = result_type::Zero();
        rv.diagonal()(0) = fix(size, source)();
        rv.diagonal()(1) = fix(size, source)();
        rv.diagonal()(2) = fix(size, source)();
        if (!VectorizedTensor::isTensorConvertibleTo1d(rv))
        {
            OGS_FATAL(
                "VectorizedTensorTest; 1d test matrix generation failed.");
        }
        return rv;
    }
};

struct random2dConvertibleEigenMatrixGenerator
{
    ac::generator<double> source;
    using result_type = Eigen::Matrix<double, 3, 3>;

    result_type operator()(std::size_t size = 0)
    {
        result_type rv = random1dConvertibleEigenMatrixGenerator{}(size);
        rv(0, 1) = fix(size, source)();
        rv(1, 0) = fix(size, source)();
        if (!VectorizedTensor::isTensorConvertibleTo2d(rv))
        {
            OGS_FATAL(
                "VectorizedTensorTest; 2d test matrix generation failed.");
        }
        return rv;
    }
};

TEST(VectorizedTensorTest, SelfTestMapping1)
{
    vectorizedTensorToMatrixMapping<1>();
    matrixToVectorizedTensorMapping<1>(
        random1dConvertibleEigenMatrixGenerator{});
}
TEST(VectorizedTensorTest, SelfTestMapping2)
{
    vectorizedTensorToMatrixMapping<2>();
    matrixToVectorizedTensorMapping<2>(
        random2dConvertibleEigenMatrixGenerator{});
}
TEST(VectorizedTensorTest, SelfTestMapping3)
{
    vectorizedTensorToMatrixMapping<3>();
    matrixToVectorizedTensorMapping<3>(
        ac::randomEigenMatrixGenerator<double, 3, 3>{});
}

TEST(VectorizedTensorTest, DynamicSizeTensorToVector)
{
    // Testing only for dimension 1, because that code part is dimension
    // independent.
    EXPECT_NO_THROW(VectorizedTensor::toVector<1>(Eigen::MatrixXd::Zero(3, 3)));
    EXPECT_THROW(VectorizedTensor::toVector<1>(Eigen::MatrixXd::Zero(3, 4)),
                 std::runtime_error);
    EXPECT_THROW(VectorizedTensor::toVector<1>(Eigen::MatrixXd::Zero(4, 3)),
                 std::runtime_error);
    EXPECT_THROW(VectorizedTensor::toVector<1>(Eigen::MatrixXd::Zero(3, 2)),
                 std::runtime_error);
    EXPECT_THROW(VectorizedTensor::toVector<1>(Eigen::MatrixXd::Zero(2, 3)),
                 std::runtime_error);
}

TEST(VectorizedTensorTest, ConvertibleAndNonconvertibleTensors)
{
    auto one_ij = [](int i, int j)
    {
        Eigen::Matrix3d m = Eigen::Matrix3d::Zero();
        m(i, j) = 1.;
        return m;
    };

    // 1d
    EXPECT_TRUE(
        VectorizedTensor::isTensorConvertibleTo1d(Eigen::MatrixXd::Zero(3, 3)));
    EXPECT_TRUE(VectorizedTensor::isTensorConvertibleTo1d(
        Eigen::MatrixXd::Identity(3, 3)));
    using namespace ranges::views;
    ranges::for_each(
        cartesian_product(iota(0, 3), iota(0, 3)),
        [&](auto const ij)
        {
            auto const [i, j] = ij;
            EXPECT_EQ(i == j,
                      VectorizedTensor::isTensorConvertibleTo1d(one_ij(i, j)));
        });

    // 2d
    EXPECT_TRUE(
        VectorizedTensor::isTensorConvertibleTo2d(Eigen::MatrixXd::Zero(3, 3)));
    EXPECT_TRUE(VectorizedTensor::isTensorConvertibleTo2d(
        Eigen::MatrixXd::Identity(3, 3)));
    {
        EXPECT_FALSE(VectorizedTensor::isTensorConvertibleTo2d(one_ij(0, 2)));
        EXPECT_FALSE(VectorizedTensor::isTensorConvertibleTo2d(one_ij(1, 2)));
        EXPECT_FALSE(VectorizedTensor::isTensorConvertibleTo2d(one_ij(2, 0)));
        EXPECT_FALSE(VectorizedTensor::isTensorConvertibleTo2d(one_ij(2, 1)));
    }
}
//
// Identities of a vectorized tensors correspond to Eigen matrix identities.
//
TEST(VectorizedTensorTest, Identity)
{
    EXPECT_EQ(VectorizedTensor::toTensor<1>(VectorizedTensor::identity<1>()),
              Eigen::Matrix3d::Identity());
    EXPECT_EQ(VectorizedTensor::toTensor<2>(VectorizedTensor::identity<2>()),
              Eigen::Matrix3d::Identity());
    EXPECT_EQ(VectorizedTensor::toTensor<3>(VectorizedTensor::identity<3>()),
              Eigen::Matrix3d::Identity());
}

//
// Determinants of a vectorized tensor and corresponding Eigen::Matrix are
// equal.
//
template <int DisplacementDim>
auto vectorizedTensorDeterminant()
{
    static constexpr int size = VectorizedTensor::size(DisplacementDim);
    ac::randomEigenMatrixGenerator<double, size, 1> generator;

    auto f = [](VectorizedTensor::Type<DisplacementDim> const& v)
    {
        return std::abs(VectorizedTensor::determinant(v) -
                        VectorizedTensor::toTensor<DisplacementDim>(v)
                            .determinant()) <=
               std::numeric_limits<double>::epsilon() *
                   std::pow(v.norm(), 3.07);
    };

    ac::check<VectorizedTensor::Type<DisplacementDim>>(
        f, 1000, ac::make_arbitrary(generator), ac::gtest_reporter());
}

TEST(VectorizedTensorTest, Determinant1)
{
    vectorizedTensorDeterminant<1>();
}
TEST(VectorizedTensorTest, Determinant2)
{
    vectorizedTensorDeterminant<2>();
}
TEST(VectorizedTensorTest, Determinant3)
{
    vectorizedTensorDeterminant<3>();
}
