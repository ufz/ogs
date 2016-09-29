/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <cmath>

#include <gtest/gtest.h>

#include <Eigen/Eigen>

#include "MeshLib/Elements/Line.h"

#include "ProcessLib/SmallDeformationWithLIE/Common/Utils.h"


namespace
{
typedef MeshLib::Line ElementType;
const unsigned e_nnodes = ElementType::n_all_nodes;

std::unique_ptr<MeshLib::Line> createLine(
    std::array<double, 3> const& a, std::array<double, 3> const& b)
{
    MeshLib::Node** nodes = new MeshLib::Node*[e_nnodes];
    nodes[0] = new MeshLib::Node(a);
    nodes[1] = new MeshLib::Node(b);
    return std::unique_ptr<MeshLib::Line>{new MeshLib::Line(nodes)};
}

std::unique_ptr<MeshLib::Line> createX()
{
    return createLine({{-1.0, 0.0, 0.0}}, {{1.0,  0.0, 0.0}});
}

std::unique_ptr<MeshLib::Line> createY()
{
    return createLine({{0.0, -1.0, 0.0}}, {{0.0,  1.0, 0.0}});
}

std::unique_ptr<MeshLib::Line> createXY()
{
    // 45degree inclined
    return createLine({{0.0, 0.0, 0.0}}, {{2./sqrt(2), 2./sqrt(2), 0.0}});
}

const double eps = std::numeric_limits<double>::epsilon();

}

TEST(LIE, rotationMatrixX)
{
    auto e(createX());
    Eigen::Vector3d nv;
    ProcessLib::SmallDeformationWithLIE::computeNormalVector(*e, nv);
    ASSERT_EQ(0., nv[0]);
    ASSERT_EQ(1., nv[1]);
    ASSERT_EQ(0., nv[2]);

    Eigen::MatrixXd R(2,2);
    ProcessLib::SmallDeformationWithLIE::computeRotationMatrix(nv, 2, R);

    ASSERT_NEAR(1., R(0,0), eps);
    ASSERT_NEAR(0., R(0,1), eps);
    ASSERT_NEAR(0., R(1,0), eps);
    ASSERT_NEAR(1., R(1,1), eps);
}

TEST(LIE, rotationMatrixY)
{
    auto e(createY());
    Eigen::Vector3d nv;
    ProcessLib::SmallDeformationWithLIE::computeNormalVector(*e, nv);
    ASSERT_EQ(-1., nv[0]);
    ASSERT_EQ(0., nv[1]);
    ASSERT_EQ(0., nv[2]);

    Eigen::MatrixXd R(2,2);
    ProcessLib::SmallDeformationWithLIE::computeRotationMatrix(nv, 2, R);

    ASSERT_NEAR(0., R(0,0), eps);
    ASSERT_NEAR(1., R(0,1), eps);
    ASSERT_NEAR(-1., R(1,0), eps);
    ASSERT_NEAR(0., R(1,1), eps);
}

TEST(LIE, rotationMatrixXY)
{
    auto e(createXY());
    Eigen::Vector3d nv;
    ProcessLib::SmallDeformationWithLIE::computeNormalVector(*e, nv);
    ASSERT_NEAR(-1./sqrt(2), nv[0], eps);
    ASSERT_NEAR(1./sqrt(2), nv[1], eps);
    ASSERT_EQ(0., nv[2]);

    Eigen::MatrixXd R(2,2);
    ProcessLib::SmallDeformationWithLIE::computeRotationMatrix(nv, 2, R);

    ASSERT_NEAR(1./sqrt(2), R(0,0), eps);
    ASSERT_NEAR(1./sqrt(2), R(0,1), eps);
    ASSERT_NEAR(-1./sqrt(2), R(1,0), eps);
    ASSERT_NEAR(1./sqrt(2), R(1,1), eps);
}
