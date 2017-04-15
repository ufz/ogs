/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <cmath>
#include <memory>

#include <gtest/gtest.h>

#include <Eigen/Eigen>

#include "MeshLib/Elements/Line.h"
#include "MeshLib/Elements/Tri.h"
#include "MeshLib/Mesh.h"

#include "ProcessLib/LIE/Common/Utils.h"


namespace
{
std::unique_ptr<MeshLib::Mesh> createTriangle(
    std::array<std::array<double, 3>, 3> const& points)
{
    auto** nodes = new MeshLib::Node*[3];
    for (int i = 0; i < points.size(); ++i)
        nodes[i] = new MeshLib::Node(points[i]);
    MeshLib::Element* e = new MeshLib::Tri(nodes);

    return std::unique_ptr<MeshLib::Mesh>(new MeshLib::Mesh(
        "", std::vector<MeshLib::Node*>{nodes[0], nodes[1], nodes[2]},
        std::vector<MeshLib::Element*>{e}));
}

std::unique_ptr<MeshLib::Mesh> createLine(
    std::array<std::array<double, 3>, 2> const& points)
{
    auto** nodes = new MeshLib::Node*[2];
    for (int i = 0; i < points.size(); ++i)
        nodes[i] = new MeshLib::Node(points[i]);
    MeshLib::Element* e = new MeshLib::Line(nodes);

    return std::unique_ptr<MeshLib::Mesh>(
        new MeshLib::Mesh("", std::vector<MeshLib::Node*>{nodes[0], nodes[1]},
                          std::vector<MeshLib::Element*>{e}));
}

const double eps = std::numeric_limits<double>::epsilon();

}

TEST(LIE, rotationMatrixXYTriangle)
{
    auto msh = createTriangle(
        {{{{0.0, 0.0, 0.0}}, {{1.0, 0.0, 0.0}}, {{1.0, 1.0, 0.0}}}});
    auto e(msh->getElement(0));
    Eigen::Vector3d nv;
    ProcessLib::LIE::computeNormalVector(*e, 3, nv);
    ASSERT_EQ(0., nv[0]);
    ASSERT_EQ(0., nv[1]);
    ASSERT_EQ(-1., nv[2]);

    Eigen::MatrixXd R(3, 3);
    ProcessLib::LIE::computeRotationMatrix(*e, nv, 3, R);

    ASSERT_NEAR(-1., R(0, 0), eps);
    ASSERT_NEAR(0., R(0, 1), eps);
    ASSERT_NEAR(0., R(0, 2), eps);
    ASSERT_NEAR(0., R(1, 0), eps);
    ASSERT_NEAR(1., R(1, 1), eps);
    ASSERT_NEAR(0., R(1, 2), eps);
    ASSERT_NEAR(0., R(2, 0), eps);
    ASSERT_NEAR(0., R(2, 1), eps);
    ASSERT_NEAR(-1., R(2, 2), eps);
}

TEST(LIE, rotationMatrixYZTriangle)
{
    auto msh = createTriangle(
        {{{{0.0, 0.0, 0.0}}, {{0.0, 1.0, 0.0}}, {{0.0, 1.0, 1.0}}}});
    auto e(msh->getElement(0));
    Eigen::Vector3d nv;
    ProcessLib::LIE::computeNormalVector(*e, 3, nv);
    ASSERT_EQ(-1., nv[0]);
    ASSERT_EQ(0., nv[1]);
    ASSERT_EQ(0., nv[2]);

    Eigen::MatrixXd R(3, 3);
    ProcessLib::LIE::computeRotationMatrix(*e, nv, 3, R);

    ASSERT_NEAR(0., R(0, 0), eps);
    ASSERT_NEAR(-1., R(0, 1), eps);
    ASSERT_NEAR(0., R(0, 2), eps);
    ASSERT_NEAR(0., R(1, 0), eps);
    ASSERT_NEAR(0., R(1, 1), eps);
    ASSERT_NEAR(1., R(1, 2), eps);
    ASSERT_NEAR(-1., R(2, 0), eps);
    ASSERT_NEAR(0., R(2, 1), eps);
    ASSERT_NEAR(0., R(2, 2), eps);
}

TEST(LIE, rotationMatrixX)
{
    auto msh = createLine({{{{-1.0, 0.0, 0.0}}, {{1.0, 0.0, 0.0}}}});
    auto e(msh->getElement(0));
    Eigen::Vector3d nv;
    ProcessLib::LIE::computeNormalVector(*e, 2, nv);
    ASSERT_EQ(0., nv[0]);
    ASSERT_EQ(1., nv[1]);
    ASSERT_EQ(0., nv[2]);

    Eigen::MatrixXd R(2,2);
    ProcessLib::LIE::computeRotationMatrix(*e, nv, 2, R);

    ASSERT_NEAR(1., R(0,0), eps);
    ASSERT_NEAR(0., R(0,1), eps);
    ASSERT_NEAR(0., R(1,0), eps);
    ASSERT_NEAR(1., R(1,1), eps);
}

TEST(LIE, rotationMatrixY)
{
    auto msh = createLine({{{{0.0, -1.0, 0.0}}, {{0.0, 1.0, 0.0}}}});
    auto e(msh->getElement(0));
    Eigen::Vector3d nv;
    ProcessLib::LIE::computeNormalVector(*e, 2, nv);
    ASSERT_EQ(-1., nv[0]);
    ASSERT_EQ(0., nv[1]);
    ASSERT_EQ(0., nv[2]);

    Eigen::MatrixXd R(2,2);
    ProcessLib::LIE::computeRotationMatrix(*e, nv, 2, R);

    ASSERT_NEAR(0., R(0,0), eps);
    ASSERT_NEAR(1., R(0,1), eps);
    ASSERT_NEAR(-1., R(1,0), eps);
    ASSERT_NEAR(0., R(1,1), eps);
}

TEST(LIE, rotationMatrixXY)
{
    // 45degree inclined
    auto msh = createLine(
        {{{{0.0, 0.0, 0.0}}, {{2. / std::sqrt(2), 2. / std::sqrt(2), 0.0}}}});
    auto e(msh->getElement(0));
    Eigen::Vector3d nv;
    ProcessLib::LIE::computeNormalVector(*e, 2, nv);
    ASSERT_NEAR(-1./std::sqrt(2), nv[0], eps);
    ASSERT_NEAR(1./std::sqrt(2), nv[1], eps);
    ASSERT_EQ(0., nv[2]);

    Eigen::MatrixXd R(2,2);
    ProcessLib::LIE::computeRotationMatrix(*e, nv, 2, R);

    ASSERT_NEAR(1./std::sqrt(2), R(0,0), eps);
    ASSERT_NEAR(1./std::sqrt(2), R(0,1), eps);
    ASSERT_NEAR(-1./std::sqrt(2), R(1,0), eps);
    ASSERT_NEAR(1./std::sqrt(2), R(1,1), eps);
}
