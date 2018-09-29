/**
 * @copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include <memory>
#include <cmath>

#include "gtest/gtest.h"

#include "MathLib/MathTools.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Elements/Hex.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshGenerators/MeshGenerator.h"
#include "MeshLib/Node.h"

using namespace MeshLib;

TEST(MeshLib, MeshGeneratorRegularHex)
{
    const double L = 10.0;
    const std::size_t n_subdivisions = 9;
    const double dL = L / static_cast<double>(n_subdivisions);
    std::unique_ptr<Mesh> msh(MeshGenerator::generateRegularHexMesh(L, n_subdivisions));

    ASSERT_EQ(std::pow(n_subdivisions, 3), msh->getNumberOfElements());
    ASSERT_EQ(std::pow(n_subdivisions+1, 3), msh->getNumberOfNodes());

    // check nodes
    const Node& node0 = *msh->getNode(0);
    const Node& node1 = *msh->getNode(1);
    const Node& node_n = *msh->getNode(msh->getNumberOfNodes()-2);
    const Node& node_n1 = *msh->getNode(msh->getNumberOfNodes()-1);
    ASSERT_DOUBLE_EQ(.0, node0[0]);
    ASSERT_DOUBLE_EQ(.0, node0[1]);
    ASSERT_DOUBLE_EQ(.0, node0[2]);
    ASSERT_DOUBLE_EQ(dL, node1[0]);
    ASSERT_DOUBLE_EQ(.0, node1[1]);
    ASSERT_DOUBLE_EQ(.0, node1[2]);
    ASSERT_DOUBLE_EQ(L-dL, node_n[0]);
    ASSERT_DOUBLE_EQ(L, node_n[1]);
    ASSERT_DOUBLE_EQ(L, node_n[2]);
    ASSERT_DOUBLE_EQ(L, node_n1[0]);
    ASSERT_DOUBLE_EQ(L, node_n1[1]);
    ASSERT_DOUBLE_EQ(L, node_n1[2]);

    // check elements
    const Element& ele0 = *msh->getElement(0);
    const Element& ele_n = *msh->getElement(msh->getNumberOfElements()-1);
    const std::size_t offset_y0 = (n_subdivisions+1);
    const std::size_t offset_z0 = (n_subdivisions+1)*(n_subdivisions+1);
    ASSERT_EQ(0u, ele0.getNodeIndex(0));
    ASSERT_EQ(1u, ele0.getNodeIndex(1));
    ASSERT_EQ(offset_y0+1, ele0.getNodeIndex(2));
    ASSERT_EQ(offset_y0, ele0.getNodeIndex(3));
    ASSERT_EQ(offset_z0, ele0.getNodeIndex(4));
    ASSERT_EQ(offset_z0+1, ele0.getNodeIndex(5));
    ASSERT_EQ(offset_z0+offset_y0+1, ele0.getNodeIndex(6));
    ASSERT_EQ(offset_z0+offset_y0, ele0.getNodeIndex(7));
    const std::size_t offset_yn0 = (n_subdivisions+1)*(n_subdivisions-1);
    const std::size_t offset_yn1 = (n_subdivisions+1)*n_subdivisions;
    const std::size_t offset_zn0 = (n_subdivisions+1)*(n_subdivisions+1)*(n_subdivisions-1);
    const std::size_t offset_zn1 = (n_subdivisions+1)*(n_subdivisions+1)*n_subdivisions;
    ASSERT_EQ(offset_zn0+offset_yn0+n_subdivisions-1, ele_n.getNodeIndex(0));
    ASSERT_EQ(offset_zn0+offset_yn0+n_subdivisions, ele_n.getNodeIndex(1));
    ASSERT_EQ(offset_zn0+offset_yn1+n_subdivisions, ele_n.getNodeIndex(2));
    ASSERT_EQ(offset_zn0+offset_yn1+n_subdivisions-1, ele_n.getNodeIndex(3));
    ASSERT_EQ(offset_zn1+offset_yn0+n_subdivisions-1, ele_n.getNodeIndex(4));
    ASSERT_EQ(offset_zn1+offset_yn0+n_subdivisions, ele_n.getNodeIndex(5));
    ASSERT_EQ(offset_zn1+offset_yn1+n_subdivisions, ele_n.getNodeIndex(6));
    ASSERT_EQ(offset_zn1+offset_yn1+n_subdivisions-1, ele_n.getNodeIndex(7));

    std::unique_ptr<Mesh> msh2 (MeshGenerator::generateRegularHexMesh(n_subdivisions, n_subdivisions, n_subdivisions, L/n_subdivisions));
    ASSERT_EQ(msh->getNumberOfNodes(), msh2->getNumberOfNodes());
    ASSERT_DOUBLE_EQ(0, MathLib::sqrDist(*(msh->getNode(msh->getNumberOfNodes()-1)), *(msh2->getNode(msh->getNumberOfNodes()-1))));

    unsigned n_x (10);
    unsigned n_y (5);
    unsigned n_z (2);
    double delta (1.2);
    std::unique_ptr<Mesh> hex_mesh (MeshGenerator::generateRegularHexMesh(n_x, n_y, n_z, delta));
    ASSERT_EQ(n_x * n_y * n_z, hex_mesh->getNumberOfElements());
    ASSERT_EQ((n_x+1) * (n_y+1) * (n_z+1), hex_mesh->getNumberOfNodes());
    const MeshLib::Node* node (hex_mesh->getNode(hex_mesh->getNumberOfNodes()-1));
    ASSERT_DOUBLE_EQ(n_x*delta, (*node)[0]);
    ASSERT_DOUBLE_EQ(n_y*delta, (*node)[1]);
    ASSERT_DOUBLE_EQ(n_z*delta, (*node)[2]);
}

TEST(MeshLib, MeshGeneratorRegularQuad)
{
    unsigned n_x (10);
    unsigned n_y (5);
    double delta (1.2);
    std::unique_ptr<Mesh> quad_mesh (MeshGenerator::generateRegularQuadMesh(n_x, n_y,delta));
    ASSERT_EQ(n_x * n_y, quad_mesh->getNumberOfElements());
    ASSERT_EQ((n_x+1) * (n_y+1), quad_mesh->getNumberOfNodes());
    const MeshLib::Node* node (quad_mesh->getNode(quad_mesh->getNumberOfNodes()-1));
    ASSERT_DOUBLE_EQ(n_x*delta, (*node)[0]);
    ASSERT_DOUBLE_EQ(n_y*delta, (*node)[1]);
    ASSERT_DOUBLE_EQ(0, (*node)[2]);

    const double L = 10.0;
    const std::size_t n_subdivisions = 9;
    std::unique_ptr<Mesh> quad_mesh2(MeshGenerator::generateRegularQuadMesh(L, n_subdivisions));
    ASSERT_EQ(n_subdivisions * n_subdivisions, quad_mesh2->getNumberOfElements());
    node = quad_mesh2->getNode(quad_mesh2->getNumberOfNodes()-1);
    ASSERT_DOUBLE_EQ(L, (*node)[0]);
    ASSERT_DOUBLE_EQ(L, (*node)[1]);
}

TEST(MeshLib, MeshGeneratorRegularPrism)
{
    double const l_x(10);
    double const l_y(10);
    double const l_z(8);
    unsigned n_x (10);
    unsigned n_y (5);
    unsigned n_z (4);
    std::unique_ptr<Mesh> mesh (MeshGenerator::generateRegularPrismMesh(l_x, l_y, l_z, n_x, n_y, n_z));
    ASSERT_EQ(2 * n_x * n_y * n_z, mesh->getNumberOfElements());
    ASSERT_EQ((n_x+1) * (n_y+1) * (n_z+1), mesh->getNumberOfNodes());

    std::size_t count (0);
    for (std::size_t k=0; k<(n_z+1); ++k)
        for (std::size_t j=0; j<(n_y+1); ++j)
            for (std::size_t i=0; i<(n_x+1); ++i)
            {
                const MeshLib::Node* node (mesh->getNode(count++));
                ASSERT_DOUBLE_EQ(static_cast<double>(i),   (*node)[0]);
                ASSERT_DOUBLE_EQ(static_cast<double>(j*2), (*node)[1]);
                ASSERT_DOUBLE_EQ(static_cast<double>(k*2), (*node)[2]);
            }
}
