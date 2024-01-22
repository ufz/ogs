/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on January 12, 2024, 3:03 PM
 */
#include <gtest/gtest.h>

#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

#include "BaseLib/ConfigTree.h"
#include "BaseLib/Logging.h"
#include "ParameterLib/ConstantParameter.h"
#include "ParameterLib/CoordinateSystem.h"
#include "ParameterLib/CreateCoordinateSystem.h"
#include "ParameterLib/SpatialPosition.h"
#include "Tests/TestTools.h"

namespace Tests
{

std::optional<ParameterLib::CoordinateSystem> createCoordinateSystem(
    std::string const& xml,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters)
{
    auto ptree = Tests::readXml(xml.c_str());
    BaseLib::ConfigTree conf(std::move(ptree), "", BaseLib::ConfigTree::onerror,
                             BaseLib::ConfigTree::onwarning);
    auto const& config =
        conf.getConfigSubtreeOptional("local_coordinate_system");

    return ParameterLib::createCoordinateSystem(config, parameters);
}
}  // namespace Tests

static std::string const parameter_name = "unit_vector";

class CoordinateSystemWithImplicitBase : public ::testing::Test
{
public:
    template <int Dimension>
    void test(std::vector<double> const& unit_vector);
};

// Test ParameterLib::CoordinateSystem::transformation_3d
void testTransformation3D(
    ParameterLib::SpatialPosition const& pos,
    ParameterLib::CoordinateSystem const& coordinate_system,
    Eigen::MatrixXd const& transform_matrix)
{
    auto const transform_matrix_3d = coordinate_system.transformation_3d(pos);

    ASSERT_NEAR((transform_matrix_3d.topLeftCorner<2, 2>() -
                 transform_matrix.topLeftCorner<2, 2>())
                    .norm(),
                0.0, std::numeric_limits<double>::epsilon());
}

template <>
void CoordinateSystemWithImplicitBase::test<2>(
    std::vector<double> const& unit_vector)
{
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> parameters;

    parameters.emplace_back(
        std::make_unique<ParameterLib::ConstantParameter<double>>(
            parameter_name, unit_vector));

    char constexpr xml[] =
        "<local_coordinate_system>"
        "   <basis_vector_0 implicit = \"true\" />"
        "   <basis_vector_1>unit_vector</basis_vector_1>"
        "</local_coordinate_system>";

    auto const coordinate_system =
        *(Tests::createCoordinateSystem(xml, parameters));

    ParameterLib::SpatialPosition pos;

    // Check CoordinateSystem::transformation<2>
    auto const transform_matrix = coordinate_system.transformation<2>(pos);

    ASSERT_NEAR((transform_matrix.transpose() * transform_matrix -
                 Eigen::Matrix<double, 2, 2>::Identity())
                    .norm(),
                0.0, std::numeric_limits<double>::epsilon());

    // Check CoordinateSystem::transformation_3d
    testTransformation3D(pos, coordinate_system, transform_matrix);
}

// Check CoordinateSystem::transformation<3>

template <>
void CoordinateSystemWithImplicitBase::test<3>(
    std::vector<double> const& unit_vector)
{
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> parameters;
    parameters.emplace_back(
        std::make_unique<ParameterLib::ConstantParameter<double>>(
            parameter_name, unit_vector));

    char constexpr xml[] =
        "<local_coordinate_system>"
        "   <basis_vector_0 implicit = \"true\" />"
        "   <basis_vector_1 implicit = \"true\" />"
        "   <basis_vector_2>unit_vector</basis_vector_2>"
        "</local_coordinate_system>";
    auto const coordinate_system =
        *(Tests::createCoordinateSystem(xml, parameters));

    ParameterLib::SpatialPosition pos;
    auto const transform_matrix = coordinate_system.transformation<3>(pos);

    double const tol = 10.0 * std::numeric_limits<double>::epsilon();
    ASSERT_NEAR((transform_matrix.transpose() * transform_matrix -
                 Eigen::Matrix<double, 3, 3>::Identity())
                    .norm(),
                0.0, tol);

    // Check CoordinateSystem::transformation_3d
    testTransformation3D(pos, coordinate_system, transform_matrix);
}

TEST_F(CoordinateSystemWithImplicitBase, test3D0)
{
    test<3>({1.0, 0.0, 0.0});
}

TEST_F(CoordinateSystemWithImplicitBase, test3D1)
{
    test<3>({0.0, 1.0, 0.0});
}

TEST_F(CoordinateSystemWithImplicitBase, test3D2)
{
    test<3>({0.0, 0.0, 1.0});
}

TEST_F(CoordinateSystemWithImplicitBase, test3D3)
{
    test<3>({-0.573576436351046, 0.0, 0.8191520442889918});
}

TEST_F(CoordinateSystemWithImplicitBase, test3D4)
{
    test<3>({0.0, -0.573576436351046, 0.8191520442889918});
}

TEST_F(CoordinateSystemWithImplicitBase, test3D5)
{
    test<3>({-0.573576436351046, 0.8191520442889918, 0.0});
}

TEST_F(CoordinateSystemWithImplicitBase, test3D6)
{
    double const alpha = 0.2;
    double const projected = std::sin(alpha);
    double const beta = -0.3;

    test<3>({std::cos(alpha), projected * std::cos(beta),
             projected * std::sin(beta)});
}

TEST_F(CoordinateSystemWithImplicitBase, test2D)
{
    double const alpha = 0.2;
    test<2>({std::cos(alpha), std::sin(alpha)});
}

TEST(CoordinateSystem, test2D)
{
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> parameters;

    double const alpha = 0.2;
    std::vector<double> e0{std::cos(alpha), std::sin(alpha)};
    parameters.emplace_back(
        std::make_unique<ParameterLib::ConstantParameter<double>>("e0", e0));
    std::vector<double> e1{-std::sin(alpha), std::cos(alpha)};
    parameters.emplace_back(
        std::make_unique<ParameterLib::ConstantParameter<double>>("e1", e1));

    char constexpr xml[] =
        "<local_coordinate_system>"
        "   <basis_vector_0>e0</basis_vector_0>"
        "   <basis_vector_1>e1</basis_vector_1>"
        "</local_coordinate_system>";
    auto const coordinate_system =
        *(Tests::createCoordinateSystem(xml, parameters));

    ParameterLib::SpatialPosition pos;
    auto const transform_matrix = coordinate_system.transformation<2>(pos);

    ASSERT_NEAR((transform_matrix.transpose() * transform_matrix -
                 Eigen::Matrix<double, 2, 2>::Identity())
                    .norm(),
                0.0, std::numeric_limits<double>::epsilon());

    // Check CoordinateSystem::transformation_3d
    testTransformation3D(pos, coordinate_system, transform_matrix);
}

TEST(CoordinateSystem, test3D)
{
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> parameters;

    std::vector<double> e0{-0.8191520442889918, 0.0, -0.573576436351046};
    parameters.emplace_back(
        std::make_unique<ParameterLib::ConstantParameter<double>>("e0", e0));
    std::vector<double> e1{0.0, -1.0, 0.0};
    parameters.emplace_back(
        std::make_unique<ParameterLib::ConstantParameter<double>>("e1", e1));
    std::vector<double> e2{-0.573576436351046, 0.0, 0.8191520442889918};
    parameters.emplace_back(
        std::make_unique<ParameterLib::ConstantParameter<double>>("e2", e2));

    char constexpr xml[] =
        "<local_coordinate_system>"
        "   <basis_vector_0>e0</basis_vector_0>"
        "   <basis_vector_1>e1</basis_vector_1>"
        "   <basis_vector_2>e2</basis_vector_2>"
        "</local_coordinate_system>";
    auto const coordinate_system =
        *(Tests::createCoordinateSystem(xml, parameters));

    ParameterLib::SpatialPosition pos;
    auto const transform_matrix = coordinate_system.transformation<3>(pos);

    ASSERT_NEAR((transform_matrix.transpose() * transform_matrix -
                 Eigen::Matrix<double, 3, 3>::Identity())
                    .norm(),
                0.0, std::numeric_limits<double>::epsilon());

    // Check CoordinateSystem::transformation_3d
    testTransformation3D(pos, coordinate_system, transform_matrix);
}
