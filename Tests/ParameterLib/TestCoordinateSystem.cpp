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

template <int Dimension>
void testException(
    std::string const& imput_xml,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    std::string_view const error_message)
{
    try
    {
        Tests::createCoordinateSystem(imput_xml.data(), parameters);
    }
    catch (std::exception& e)
    {
        EXPECT_EQ(e.what(), error_message)
            << "Expected error message:\n"
            << e.what() << "\nGiven error message:\n"
            << error_message;
    }
}

static std::string const error_message_implicit_explicit_2D =
    "The case of implicit \"basis_vector_0\" and explicit "
    "\"basis_vector_1\" is for a 2D coordinate system. The parameter "
    "for \"basis_vector_1\", e1, must have two components but it has "
    "3. In addition, \"basis_vector_2\" should not exist in this case.";

static std::string const error_message_explicit_implicit =
    "Since basis_vector_0 is explicitly defined, basis_vector_1 must "
    "be explicit as well.";

static std::string const error_message_implicit_base2 =
    "basis_vector_2 must be explicit.";

TEST(CoordinateSystem, testWrongBase0ComponentNumberException)
{
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> parameters;

    std::vector<double> e0{1.0, 0.0, 0.0, 0.0};
    parameters.emplace_back(
        std::make_unique<ParameterLib::ConstantParameter<double>>("e0", e0));

    std::string const error_message =
        "Basis vector parameter 'e0' must have two or three components, but it "
        "has 4.";

    std::string const xml =
        "<local_coordinate_system>"
        "   <basis_vector_0>e0</basis_vector_0>"
        "</local_coordinate_system>";

    testException<3>(xml, parameters, error_message);
}

// Test: basis_vector_2 must not be presented in  the 2D coordinate system

TEST(CoordinateSystem, test2DBase2Exception)
{
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> parameters;

    std::vector<double> e0{1.0, 0.0};
    parameters.emplace_back(
        std::make_unique<ParameterLib::ConstantParameter<double>>("e0", e0));
    std::vector<double> e1{0.0, 1.0};
    parameters.emplace_back(
        std::make_unique<ParameterLib::ConstantParameter<double>>("e1", e1));

    std::string const error_message =
        "The tag \"basis_vector_2\" is not needed for a 2D local coordinate "
        "system.";

    {
        std::string const xml =
            "<local_coordinate_system>"
            "   <basis_vector_0>e0</basis_vector_0>"
            "   <basis_vector_1>e1</basis_vector_1>"
            "   <basis_vector_2>e1</basis_vector_2>"
            "</local_coordinate_system>";
        testException<2>(xml, parameters, error_message);
    }
    {
        std::string const xml =
            "<local_coordinate_system>"
            "   <basis_vector_0 implicit = \"true\" />"
            "   <basis_vector_1>e1</basis_vector_1>"
            "   <basis_vector_2 implicit = \"true\" />"
            "</local_coordinate_system>";
        testException<2>(xml, parameters, error_message);
    }
}

TEST(CoordinateSystem, test2DWrongComponentNumberException)
{
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> parameters;

    std::vector<double> e0{1.0, 0.0};
    parameters.emplace_back(
        std::make_unique<ParameterLib::ConstantParameter<double>>("e0", e0));
    std::vector<double> e1{0.0, 1.0, 0.0};
    parameters.emplace_back(
        std::make_unique<ParameterLib::ConstantParameter<double>>("e1", e1));

    {
        std::string const xml =
            "<local_coordinate_system>"
            "   <basis_vector_0>e0</basis_vector_0>"
            "   <basis_vector_1>e1</basis_vector_1>"
            "</local_coordinate_system>";
        std::string const error_message =
            "The read parameter `e1' has the wrong number of components"
            " (3 instead of 2).";

        testException<2>(xml, parameters, error_message);
    }

    {  // With implicit base

        std::string const xml =
            "<local_coordinate_system>"
            "   <basis_vector_0 implicit = \"true\" />"
            "   <basis_vector_1>e1</basis_vector_1>"
            "</local_coordinate_system>";

        testException<2>(xml, parameters, error_message_implicit_explicit_2D);
    }
}

TEST(CoordinateSystem, test3DWrongComponentNumberException)
{
    std::vector<double> e0{1.0, 0.0, 0.0};
    {  // Wrong component number in base 1
        std::vector<std::unique_ptr<ParameterLib::ParameterBase>> parameters;

        parameters.emplace_back(
            std::make_unique<ParameterLib::ConstantParameter<double>>("e0",
                                                                      e0));
        std::vector<double> e1{0.0, 1.0};
        parameters.emplace_back(
            std::make_unique<ParameterLib::ConstantParameter<double>>("e1",
                                                                      e1));
        std::vector<double> e2{0.0, 1.0, 1.0};
        parameters.emplace_back(
            std::make_unique<ParameterLib::ConstantParameter<double>>("e2",
                                                                      e2));

        std::string const xml =
            "<local_coordinate_system>"
            "   <basis_vector_0>e0</basis_vector_0>"
            "   <basis_vector_1>e1</basis_vector_1>"
            "   <basis_vector_2>e2</basis_vector_2>"
            "</local_coordinate_system>";
        std::string const error_message =
            "The read parameter `e1' has the wrong number of components"
            " (2 instead of 3).";

        testException<3>(xml, parameters, error_message);
    }

    {  // Wrong component number in base 2
        std::vector<std::unique_ptr<ParameterLib::ParameterBase>> parameters;

        parameters.emplace_back(
            std::make_unique<ParameterLib::ConstantParameter<double>>("e0",
                                                                      e0));

        std::vector<double> e1{0.0, 1.0, 0.0};
        parameters.emplace_back(
            std::make_unique<ParameterLib::ConstantParameter<double>>("e1",
                                                                      e1));
        std::vector<double> e2{0.0, 1.0};
        parameters.emplace_back(
            std::make_unique<ParameterLib::ConstantParameter<double>>("e2",
                                                                      e2));

        std::string const error_message =
            "The read parameter `e2' for tag basis_vector_2 has the wrong "
            "number of components (2 instead of 3).";

        {
            std::string const xml =
                "<local_coordinate_system>"
                "   <basis_vector_0>e0</basis_vector_0>"
                "   <basis_vector_1>e1</basis_vector_1>"
                "   <basis_vector_2>e2</basis_vector_2>"
                "</local_coordinate_system>";

            testException<3>(xml, parameters, error_message);
        }

        {  // with implicit bases

            std::string const xml =
                "<local_coordinate_system>"
                "   <basis_vector_0 implicit = \"true\" />"
                "   <basis_vector_1 implicit = \"true\" />"
                "   <basis_vector_2>e2</basis_vector_2>"
                "</local_coordinate_system>";

            testException<3>(xml, parameters, error_message);
        }
    }
}

TEST(CoordinateSystem, test2DWrongImplicitBaseException)
{
    {  // Two bases are implicit
        std::vector<std::unique_ptr<ParameterLib::ParameterBase>> parameters;

        std::string const xml =
            "<local_coordinate_system>"
            "   <basis_vector_0 implicit = \"true\" />"
            "   <basis_vector_1 implicit = \"true\" />"
            "</local_coordinate_system>";
        std::string const error_message =
            "Both \"basis_vector_0\" and \"basis_vector_1\" are implicit but "
            "\"basis_vector_2\" does not exist. If 2D coordinate system is "
            "considered, please change \"basis_vector_1\" to explicit. If 3D "
            "coordinate system is considered, please add \"basis_vector_2\".";

        testException<2>(xml, parameters, error_message);
    }

    {  // base0 is explicit but base1 is implicit

        std::vector<std::unique_ptr<ParameterLib::ParameterBase>> parameters;

        std::vector<double> e0{1.0, 0.0};
        parameters.emplace_back(
            std::make_unique<ParameterLib::ConstantParameter<double>>("e0",
                                                                      e0));

        std::string const xml =
            "<local_coordinate_system>"
            "   <basis_vector_0>e0</basis_vector_0>"
            "   <basis_vector_1 implicit = \"true\" />"
            "</local_coordinate_system>";

        testException<2>(xml, parameters, error_message_explicit_implicit);
    }
}

TEST(CoordinateSystem, test3DWrongImplicitBaseException)
{
    // Base attribute: i for implicit and e for explicit
    // case: e e i
    {  // base0 and base1 are explicit but base2 is implicit

        std::vector<std::unique_ptr<ParameterLib::ParameterBase>> parameters;

        std::vector<double> e0{1.0, 0.0, 0.0};
        parameters.emplace_back(
            std::make_unique<ParameterLib::ConstantParameter<double>>("e0",
                                                                      e0));
        std::vector<double> e1{0.0, 1.0, 0.0};
        parameters.emplace_back(
            std::make_unique<ParameterLib::ConstantParameter<double>>("e1",
                                                                      e1));

        std::string const xml =
            "<local_coordinate_system>"
            "   <basis_vector_0>e0</basis_vector_0>"
            "   <basis_vector_1>e0</basis_vector_1>"
            "   <basis_vector_2 implicit = \"true\" />"
            "</local_coordinate_system>";

        testException<3>(xml, parameters, error_message_implicit_base2);
    }

    // case: e i e
    {
        std::vector<std::unique_ptr<ParameterLib::ParameterBase>> parameters;

        std::vector<double> e0{1.0, 0.0, 0.0};
        parameters.emplace_back(
            std::make_unique<ParameterLib::ConstantParameter<double>>("e0",
                                                                      e0));
        std::vector<double> e2{0.0, 0.0, 1.0};
        parameters.emplace_back(
            std::make_unique<ParameterLib::ConstantParameter<double>>("e2",
                                                                      e2));
        std::string const xml =
            "<local_coordinate_system>"
            "   <basis_vector_0>e0</basis_vector_0>"
            "   <basis_vector_1 implicit = \"true\" />"
            "   <basis_vector_2>e2</basis_vector_2>"
            "</local_coordinate_system>";

        testException<3>(xml, parameters, error_message_explicit_implicit);
    }

    // case: e i i
    {
        std::vector<std::unique_ptr<ParameterLib::ParameterBase>> parameters;

        std::vector<double> e0{1.0, 0.0, 0.0};
        parameters.emplace_back(
            std::make_unique<ParameterLib::ConstantParameter<double>>("e0",
                                                                      e0));
        std::string const xml =
            "<local_coordinate_system>"
            "   <basis_vector_0>e0</basis_vector_0>"
            "   <basis_vector_1 implicit = \"true\" />"
            "   <basis_vector_2 implicit = \"true\" />"
            "</local_coordinate_system>";

        testException<3>(xml, parameters, error_message_explicit_implicit);
    }

    // case: i i i
    {  // all bases are implicit
        std::vector<std::unique_ptr<ParameterLib::ParameterBase>> parameters;

        std::string const xml =
            "<local_coordinate_system>"
            "   <basis_vector_0 implicit = \"true\" />"
            "   <basis_vector_1 implicit = \"true\" />"
            "   <basis_vector_2 implicit = \"true\" />"
            "</local_coordinate_system>";

        testException<3>(xml, parameters, error_message_implicit_base2);
    }

    // case: i e i
    {  // all bases are implicit
        std::vector<std::unique_ptr<ParameterLib::ParameterBase>> parameters;

        std::vector<double> e1{0.0, 1.0, 0.0};
        parameters.emplace_back(
            std::make_unique<ParameterLib::ConstantParameter<double>>("e1",
                                                                      e1));
        std::string const xml =
            "<local_coordinate_system>"
            "   <basis_vector_0 implicit = \"true\" />"
            "   <basis_vector_1>e1</basis_vector_1>"
            "   <basis_vector_2 implicit = \"true\" />"
            "</local_coordinate_system>";

        testException<3>(xml, parameters, error_message_implicit_explicit_2D);
    }

    // case: i e e
    {
        std::vector<std::unique_ptr<ParameterLib::ParameterBase>> parameters;

        std::vector<double> e1{0.0, 1.0, 0.0};
        parameters.emplace_back(
            std::make_unique<ParameterLib::ConstantParameter<double>>("e1",
                                                                      e1));
        std::vector<double> e2{0.0, 0.0, 1.0};
        parameters.emplace_back(
            std::make_unique<ParameterLib::ConstantParameter<double>>("e2",
                                                                      e2));
        std::string const xml =
            "<local_coordinate_system>"
            "   <basis_vector_0 implicit = \"true\" />"
            "   <basis_vector_1>e1</basis_vector_1>"
            "   <basis_vector_2>e2</basis_vector_2>"
            "</local_coordinate_system>";

        testException<3>(xml, parameters, error_message_implicit_explicit_2D);

        // For this case, if there is no exception for base1, i.e. with two
        // components, an exception will be caused due to the existence of
        // basis_vector_2. This exception has already been tested in
        // CoordinateSystem.test2DBase2Exception.
    }
}
