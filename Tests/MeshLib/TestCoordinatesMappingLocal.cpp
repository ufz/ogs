/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include <gtest/gtest.h>

#include <limits>
#include <algorithm>
#include <vector>
#include <cmath>
#include <memory>

#include <Eigen/Eigen>

#include "GeoLib/AnalyticalGeometry.h"
#include "MathLib/LinAlg/Dense/DenseMatrix.h"

#include "MeshLib/CoordinateSystem.h"
#include "MeshLib/Node.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Elements/Line.h"
#include "MeshLib/Elements/Quad.h"
#include "MeshLib/ElementCoordinatesMappingLocal.h"

#include "Tests/TestTools.h"


namespace
{

namespace TestLine2
{
using ElementType = MeshLib::Line;
const unsigned e_nnodes = ElementType::n_all_nodes;

std::unique_ptr<MeshLib::Line> createLine(std::array<double, 3> const& a,
                                          std::array<double, 3> const& b)
{
    auto** nodes = new MeshLib::Node*[e_nnodes];
    nodes[0] = new MeshLib::Node(a);
    nodes[1] = new MeshLib::Node(b);
    return std::unique_ptr<MeshLib::Line>{new MeshLib::Line(nodes)};
    }

    std::unique_ptr<MeshLib::Line> createY()
    {
        return createLine({{0.0, -1.0, 0.0}}, {{0.0,  1.0, 0.0}});
    }

    std::unique_ptr<MeshLib::Line> createZ()
    {
        return createLine({{0.0, 0.0, -1.0}}, {{0.0, 0.0,  1.0}});
    }

    std::unique_ptr<MeshLib::Line> createXY()
    {
        // 45degree inclined
        return createLine({{0.0, 0.0, 0.0}}, {{2./sqrt(2), 2./sqrt(2), 0.0}});
    }

    std::unique_ptr<MeshLib::Line> createXYZ()
    {
        return createLine({{0.0, 0.0, 0.0}}, {{2./sqrt(3), 2./sqrt(3), 2./sqrt(3)}});
    }

};

namespace TestQuad4
{
    // Element information
using ElementType = MeshLib::Quad;
const unsigned e_nnodes = ElementType::n_all_nodes;

std::unique_ptr<MeshLib::Quad> createQuad(std::array<double, 3> const& a,
                                          std::array<double, 3> const& b,
                                          std::array<double, 3> const& c,
                                          std::array<double, 3> const& d)
{
    auto** nodes = new MeshLib::Node*[e_nnodes];
    nodes[0] = new MeshLib::Node(a);
    nodes[1] = new MeshLib::Node(b);
    nodes[2] = new MeshLib::Node(c);
    nodes[3] = new MeshLib::Node(d);
    return std::unique_ptr<MeshLib::Quad>{new MeshLib::Quad(nodes)};
    }

    // 2.5D case: inclined
    std::unique_ptr<MeshLib::Quad> createXYZ()
    {
        // rotate 45 degree around x axis
        return createQuad(
            {{ 1.0,  0.7071067811865475,  0.7071067811865475}},
            {{-1.0,  0.7071067811865475,  0.7071067811865475}},
            {{-1.0, -0.7071067811865475, -0.7071067811865475}},
            {{ 1.0, -0.7071067811865475, -0.7071067811865475}});
    }

    // 2.5D case: inclined
    std::unique_ptr<MeshLib::Quad> createXZ()
    {
        return createQuad(
            {{ 1.0, 0.0,  1.0}},
            {{-1.0, 0.0,  1.0}},
            {{-1.0, 0.0, -1.0}},
            {{ 1.0, 0.0, -1.0}});
    }

    // 2.5D case: inclined
    std::unique_ptr<MeshLib::Quad> createYZ()
    {
        return createQuad(
            {{0.0,  1.0,  1.0}},
            {{0.0, -1.0,  1.0}},
            {{0.0, -1.0, -1.0}},
            {{0.0,  1.0, -1.0}});
    }
};

#if 0
// keep this function for debugging
void debugOutput(MeshLib::Element *ele, MeshLib::ElementCoordinatesMappingLocal &mapping)
{
    std::cout.precision(12);
    std::cout << "original" << std::endl;
    for (unsigned i=0; i<ele->getNumberOfNodes(); i++)
        std::cout << *ele->getNode(i) << std::endl;
    std::cout << "local coords=" << std::endl;
    for (unsigned i=0; i<ele->getNumberOfNodes(); i++)
        std::cout << *mapping.getMappedCoordinates(i) << std::endl;
    std::cout << "R=\n" << mapping.getRotationMatrixToGlobal() << std::endl;
    auto matR(mapping.getRotationMatrixToGlobal());
    std::cout << "global coords=" << std::endl;
    for (unsigned i=0; i<ele->getNumberOfNodes(); i++) {
        double* raw = const_cast<double*>(&(*mapping.getMappedCoordinates(i))[0]);
        Eigen::Map<Eigen::Vector3d> v(raw);
        std::cout << (matR*v).transpose() << std::endl;
    }
}
#endif

// check if using the rotation matrix results in the original coordinates
#define CHECK_COORDS(ele, mapping)\
    for (unsigned ii=0; ii<(ele)->getNumberOfNodes(); ii++) {\
        MathLib::Point3d global(matR*(mapping).getMappedCoordinates(ii));\
        const double eps(std::numeric_limits<double>::epsilon());\
        ASSERT_ARRAY_NEAR(&(*(ele)->getNode(ii))[0], global.getCoords(), 3u, eps);\
    }

} //namespace

TEST(MeshLib, CoordinatesMappingLocalLowerDimLineY)
{
    auto ele = TestLine2::createY();
    MeshLib::ElementCoordinatesMappingLocal mapping(
        *ele, MeshLib::CoordinateSystem(MeshLib::CoordinateSystemType::Y)
                  .getDimension());
    auto matR(mapping.getRotationMatrixToGlobal());
    //debugOutput(ele, mapping);

    double exp_R[3*3] = {0, -1, 0,
                         1,  0, 0,
                         0,  0, 1};
    const double eps(std::numeric_limits<double>::epsilon());
    ASSERT_ARRAY_NEAR(exp_R, matR.data(), matR.size(), eps);
    CHECK_COORDS(ele,mapping);

    for (std::size_t n = 0; n < ele->getNumberOfNodes(); ++n)
        delete ele->getNode(n);
}

TEST(MeshLib, CoordinatesMappingLocalLowerDimLineZ)
{
    auto ele = TestLine2::createZ();
    MeshLib::ElementCoordinatesMappingLocal mapping(
        *ele,
        MeshLib::CoordinateSystem(MeshLib::CoordinateSystemType::Z)
            .getDimension());
    auto matR(mapping.getRotationMatrixToGlobal());
    //debugOutput(ele, mapping);

    double exp_R[3*3] = {0, 0, -1, 0, 1, 0, 1, 0, 0};
    const double eps(std::numeric_limits<double>::epsilon());
    ASSERT_ARRAY_NEAR(exp_R, matR.data(), matR.size(), eps);
    CHECK_COORDS(ele,mapping);

    for (std::size_t n = 0; n < ele->getNumberOfNodes(); ++n)
        delete ele->getNode(n);
}

TEST(MeshLib, CoordinatesMappingLocalLowerDimLineXY)
{
    auto ele = TestLine2::createXY();
    MeshLib::ElementCoordinatesMappingLocal mapping(
        *ele, MeshLib::CoordinateSystem(*ele).getDimension());
    auto matR(mapping.getRotationMatrixToGlobal());
    //debugOutput(ele, mapping);

    double exp_R[3*3] = {0.70710678118654757, -0.70710678118654757, 0,
                         0.70710678118654757,  0.70710678118654757, 0,
                         0,                    0,                   1};
    const double eps(std::numeric_limits<double>::epsilon());
    ASSERT_ARRAY_NEAR(exp_R, matR.data(), matR.size(), eps);
    CHECK_COORDS(ele,mapping);

    for (std::size_t n = 0; n < ele->getNumberOfNodes(); ++n)
        delete ele->getNode(n);
}

TEST(MeshLib, CoordinatesMappingLocalLowerDimLineXYZ)
{
    auto ele = TestLine2::createXYZ();
    MeshLib::ElementCoordinatesMappingLocal mapping(
        *ele, MeshLib::CoordinateSystem(*ele).getDimension());
    auto matR(mapping.getRotationMatrixToGlobal());
    //debugOutput(ele, mapping);

    double exp_R[3*3] = {0.57735026918962584, -0.81649658092772626,  0,
                         0.57735026918962584,  0.40824829046386313, -0.70710678118654757,
                         0.57735026918962584,  0.40824829046386313,  0.70710678118654757};
    const double eps(std::numeric_limits<double>::epsilon());
    ASSERT_ARRAY_NEAR(exp_R, matR.data(), matR.size(), eps);
    CHECK_COORDS(ele,mapping);

    for (std::size_t n = 0; n < ele->getNumberOfNodes(); ++n)
        delete ele->getNode(n);
}

TEST(MeshLib, CoordinatesMappingLocalLowerDimQuadXZ)
{
    auto ele = TestQuad4::createXZ();
    MeshLib::ElementCoordinatesMappingLocal mapping(
        *ele, MeshLib::CoordinateSystem(*ele).getDimension());
    auto matR(mapping.getRotationMatrixToGlobal());
    //debugOutput(ele, mapping);

    // results when using GeoLib::ComputeRotationMatrixToXY()
    double exp_R[3*3] = {  1, 0,  0,
                           0, 0, -1,
                           0, 1,  0};

    const double eps(std::numeric_limits<double>::epsilon());
    ASSERT_ARRAY_NEAR(exp_R, matR.data(), matR.size(), eps);
    CHECK_COORDS(ele,mapping);

    for (std::size_t n = 0; n < ele->getNumberOfNodes(); ++n)
        delete ele->getNode(n);
}

TEST(MeshLib, CoordinatesMappingLocalLowerDimQuadYZ)
{
    auto ele = TestQuad4::createYZ();
    MeshLib::ElementCoordinatesMappingLocal mapping(
        *ele, MeshLib::CoordinateSystem(*ele).getDimension());
    auto matR(mapping.getRotationMatrixToGlobal());
    //debugOutput(ele, mapping);

    // results when using GeoLib::ComputeRotationMatrixToXY()
    double exp_R[3*3] = { 0, 0, 1,
                          0, 1, 0,
                         -1, 0, 0};

    const double eps(std::numeric_limits<double>::epsilon());
    ASSERT_ARRAY_NEAR(exp_R, matR.data(), matR.size(), eps);
    CHECK_COORDS(ele,mapping);

    for (std::size_t n = 0; n < ele->getNumberOfNodes(); ++n)
        delete ele->getNode(n);
}

TEST(MeshLib, CoordinatesMappingLocalLowerDimQuadXYZ)
{
    auto ele = TestQuad4::createXYZ();
    MeshLib::ElementCoordinatesMappingLocal mapping(
        *ele, MeshLib::CoordinateSystem(*ele).getDimension());
    auto matR(mapping.getRotationMatrixToGlobal());
    //debugOutput(ele, mapping);

    // results when using GeoLib::ComputeRotationMatrixToXY()
    double exp_R[3*3] = {  1, 0, 0,
                           0, 0.70710678118654757, -0.70710678118654757,
                           0, 0.70710678118654757,  0.70710678118654757};

    const double eps(std::numeric_limits<double>::epsilon());
    ASSERT_ARRAY_NEAR(exp_R, matR.data(), matR.size(), eps);
    CHECK_COORDS(ele,mapping);

    for (std::size_t n = 0; n < ele->getNumberOfNodes(); ++n)
        delete ele->getNode(n);
}

