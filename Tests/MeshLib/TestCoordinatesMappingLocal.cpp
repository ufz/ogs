/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include <gtest/gtest.h>

#include <limits>
#include <algorithm>
#include <vector>
#include <cmath>

#ifdef OGS_USE_EIGEN
#include <Eigen/Eigen>
#endif

#include "GeoLib/AnalyticalGeometry.h"
#include "MathLib/LinAlg/Dense/DenseMatrix.h"

#include "MeshLib/Node.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Elements/Line.h"
#include "MeshLib/Elements/Quad.h"
#include "MeshLib/ElementCoordinatesMappingLocal.h"

#include "../TestTools.h"


namespace
{

class TestLine2
{
public:
    typedef MeshLib::Line ElementType;
    static const unsigned e_nnodes = ElementType::n_all_nodes;

    static MeshLib::Line* createY()
    {
        MeshLib::Node** nodes = new MeshLib::Node*[e_nnodes];
        nodes[0] = new MeshLib::Node(0.0, -1.0, 0.0);
        nodes[1] = new MeshLib::Node(0.0,  1.0, 0.0);
        return new MeshLib::Line(nodes);
    }

    static MeshLib::Line* createZ()
    {
        MeshLib::Node** nodes = new MeshLib::Node*[e_nnodes];
        nodes[0] = new MeshLib::Node(0.0, 0.0, -1.0);
        nodes[1] = new MeshLib::Node(0.0, 0.0,  1.0);
        return new MeshLib::Line(nodes);
    }

    static MeshLib::Line* createXY()
    {
        // 45degree inclined
        MeshLib::Node** nodes = new MeshLib::Node*[e_nnodes];
        nodes[0] = new MeshLib::Node(0.0, 0.0, 0.0);
        nodes[1] = new MeshLib::Node(2./sqrt(2), 2./sqrt(2), 0.0);
        return new MeshLib::Line(nodes);
    }

    static MeshLib::Line* createXYZ()
    {
        MeshLib::Node** nodes = new MeshLib::Node*[e_nnodes];
        nodes[0] = new MeshLib::Node(0.0, 0.0, 0.0);
        nodes[1] = new MeshLib::Node(2./sqrt(3), 2./sqrt(3), 2./sqrt(3));
        return new MeshLib::Line(nodes);
    }

};

class TestQuad4
{
 public:
    // Element information
    typedef MeshLib::Quad ElementType;
    static const unsigned e_nnodes = ElementType::n_all_nodes;

    // 2.5D case: inclined
    static MeshLib::Quad* createXYZ()
    {
        MeshLib::Node** nodes = new MeshLib::Node*[e_nnodes];
        // rotate 45 degree around x axis
        nodes[0] = new MeshLib::Node( 1.0,  0.7071067811865475,  0.7071067811865475);
        nodes[1] = new MeshLib::Node(-1.0,  0.7071067811865475,  0.7071067811865475);
        nodes[2] = new MeshLib::Node(-1.0, -0.7071067811865475, -0.7071067811865475);
        nodes[3] = new MeshLib::Node( 1.0, -0.7071067811865475, -0.7071067811865475);
        return new MeshLib::Quad(nodes);
    }

    // 2.5D case: inclined
    static MeshLib::Quad* createXZ()
    {
        MeshLib::Node** nodes = new MeshLib::Node*[e_nnodes];
        nodes[0] = new MeshLib::Node( 1.0, 0.0,  1.0);
        nodes[1] = new MeshLib::Node(-1.0, 0.0,  1.0);
        nodes[2] = new MeshLib::Node(-1.0, 0.0, -1.0);
        nodes[3] = new MeshLib::Node( 1.0, 0.0, -1.0);
        return new MeshLib::Quad(nodes);
    }

    // 2.5D case: inclined
    static MeshLib::Quad* createYZ()
    {
        MeshLib::Node** nodes = new MeshLib::Node*[e_nnodes];
        nodes[0] = new MeshLib::Node(0.0,  1.0,  1.0);
        nodes[1] = new MeshLib::Node(0.0, -1.0,  1.0);
        nodes[2] = new MeshLib::Node(0.0, -1.0, -1.0);
        nodes[3] = new MeshLib::Node(0.0,  1.0, -1.0);
        return new MeshLib::Quad(nodes);
    }

};

#if 0
// keep this function for debugging
void debugOutput(MeshLib::Element *ele, MeshLib::ElementCoordinatesMappingLocal &mapping)
{
    std::cout.precision(12);
    std::cout << "original" << std::endl;
    for (unsigned i=0; i<ele->getNNodes(); i++)
        std::cout << *ele->getNode(i) << std::endl;
    std::cout << "local coords=" << std::endl;
    for (unsigned i=0; i<ele->getNNodes(); i++)
        std::cout << *mapping.getMappedCoordinates(i) << std::endl;
    std::cout << "R=\n" << mapping.getRotationMatrixToGlobal() << std::endl;
    auto matR(mapping.getRotationMatrixToGlobal());
    std::cout << "global coords=" << std::endl;
    for (unsigned i=0; i<ele->getNNodes(); i++) {
        double* raw = const_cast<double*>(&(*mapping.getMappedCoordinates(i))[0]);
        Eigen::Map<Eigen::Vector3d> v(raw);
        std::cout << (matR*v).transpose() << std::endl;
    }
}
#endif

// check if using the rotation matrix results in the original coordinates
#define CHECK_COORDS(ele, mapping)\
    for (unsigned ii=0; ii<ele->getNNodes(); ii++) {\
        Eigen::Map<Eigen::Vector3d> local(const_cast<double*>(&(*mapping.getMappedCoordinates(ii))[0]));\
        Eigen::Vector3d global(matR*local);\
        const double eps(std::numeric_limits<double>::epsilon());\
        ASSERT_ARRAY_NEAR(&(*ele->getNode(ii))[0], global.data(), 3u, eps);\
    }

} //namespace

TEST(MeshLib, CoordinatesMappingLocalLowerDimLineY)
{
    auto ele = TestLine2::createY();
    MeshLib::ElementCoordinatesMappingLocal mapping(*ele, MeshLib::CoordinateSystem(MeshLib::CoordinateSystemType::Y));
    auto matR(mapping.getRotationMatrixToGlobal());
    //debugOutput(ele, mapping);

    double exp_R[3*3] = {0, -1, 0,
                         1,  0, 0,
                         0,  0, 1};
    const double eps(std::numeric_limits<double>::epsilon());
    ASSERT_ARRAY_NEAR(exp_R, matR.data(), matR.size(), eps);
    CHECK_COORDS(ele,mapping);
}

TEST(MeshLib, CoordinatesMappingLocalLowerDimLineZ)
{
    auto ele = TestLine2::createZ();
    MeshLib::ElementCoordinatesMappingLocal mapping(*ele,
            MeshLib::CoordinateSystem(MeshLib::CoordinateSystemType::Z));
    auto matR(mapping.getRotationMatrixToGlobal());
    //debugOutput(ele, mapping);

    double exp_R[3*3] = {0, 0, -1, 0, 1, 0, 1, 0, 0};
    const double eps(std::numeric_limits<double>::epsilon());
    ASSERT_ARRAY_NEAR(exp_R, matR.data(), matR.size(), eps);
    CHECK_COORDS(ele,mapping);
}

TEST(MeshLib, CoordinatesMappingLocalLowerDimLineXY)
{
    auto ele = TestLine2::createXY();
    MeshLib::ElementCoordinatesMappingLocal mapping(*ele, MeshLib::CoordinateSystem(*ele));
    auto matR(mapping.getRotationMatrixToGlobal());
    //debugOutput(ele, mapping);

    double exp_R[3*3] = {0.70710678118654757, -0.70710678118654757, 0,
                         0.70710678118654757,  0.70710678118654757, 0,
                         0,                    0,                   1};
    const double eps(std::numeric_limits<double>::epsilon());
    ASSERT_ARRAY_NEAR(exp_R, matR.data(), matR.size(), eps);
    CHECK_COORDS(ele,mapping);
}

TEST(MeshLib, CoordinatesMappingLocalLowerDimLineXYZ)
{
    auto ele = TestLine2::createXYZ();
    MeshLib::ElementCoordinatesMappingLocal mapping(*ele, MeshLib::CoordinateSystem(*ele));
    auto matR(mapping.getRotationMatrixToGlobal());
    //debugOutput(ele, mapping);

    double exp_R[3*3] = {0.57735026918962584, -0.81649658092772626,  0,
                         0.57735026918962584,  0.40824829046386313, -0.70710678118654757,
                         0.57735026918962584,  0.40824829046386313,  0.70710678118654757};
    const double eps(std::numeric_limits<double>::epsilon());
    ASSERT_ARRAY_NEAR(exp_R, matR.data(), matR.size(), eps);
    CHECK_COORDS(ele,mapping);
}

TEST(MeshLib, CoordinatesMappingLocalLowerDimQuadXZ)
{
    auto ele = TestQuad4::createXZ();
    MeshLib::ElementCoordinatesMappingLocal mapping(*ele, MeshLib::CoordinateSystem(*ele));
    auto matR(mapping.getRotationMatrixToGlobal());
    //debugOutput(ele, mapping);

    // results when using GeoLib::ComputeRotationMatrixToXY()
    double exp_R[3*3] = {  1, 0,  0,
                           0, 0, -1,
                           0, 1,  0};

    const double eps(std::numeric_limits<double>::epsilon());
    ASSERT_ARRAY_NEAR(exp_R, matR.data(), matR.size(), eps);
    CHECK_COORDS(ele,mapping);
}

TEST(MeshLib, CoordinatesMappingLocalLowerDimQuadYZ)
{
    auto ele = TestQuad4::createYZ();
    MeshLib::ElementCoordinatesMappingLocal mapping(*ele, MeshLib::CoordinateSystem(*ele));
    auto matR(mapping.getRotationMatrixToGlobal());
    //debugOutput(ele, mapping);

    // results when using GeoLib::ComputeRotationMatrixToXY()
    double exp_R[3*3] = { 0, 0, 1,
                          0, 1, 0,
                         -1, 0, 0};

    const double eps(std::numeric_limits<double>::epsilon());
    ASSERT_ARRAY_NEAR(exp_R, matR.data(), matR.size(), eps);
    CHECK_COORDS(ele,mapping);
}

TEST(MeshLib, CoordinatesMappingLocalLowerDimQuadXYZ)
{
    auto ele = TestQuad4::createXYZ();
    MeshLib::ElementCoordinatesMappingLocal mapping(*ele, MeshLib::CoordinateSystem(*ele));
    auto matR(mapping.getRotationMatrixToGlobal());
    //debugOutput(ele, mapping);

    // results when using GeoLib::ComputeRotationMatrixToXY()
    double exp_R[3*3] = {  1, 0, 0,
                           0, 0.70710678118654757, -0.70710678118654757,
                           0, 0.70710678118654757,  0.70710678118654757};

    const double eps(std::numeric_limits<double>::epsilon());
    ASSERT_ARRAY_NEAR(exp_R, matR.data(), matR.size(), eps);
    CHECK_COORDS(ele,mapping);
}

