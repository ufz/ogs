/**
 * @file TestMeshRevision.cpp
 * @author Karsten Rink
 * @date 2013-03-04
 * @brief Tests for MeshRevision class
 *
 * @copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include "gtest/gtest.h"

#include "MeshLib/Elements/Element.h"
#include "MeshLib/Elements/Hex.h"
#include "MeshLib/Elements/Prism.h"
#include "MeshLib/Elements/Pyramid.h"
#include "MeshLib/Elements/Quad.h"
#include "MeshLib/Elements/Tet.h"
#include "MeshLib/Elements/Tri.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshEditing/MeshRevision.h"
#include "MeshLib/Node.h"


TEST(MeshEditing, Tri)
{
    std::vector<MeshLib::Node*> nodes;
    nodes.push_back(new MeshLib::Node(1,0,0));
    nodes.push_back(new MeshLib::Node(0,0,0));
    nodes.push_back(new MeshLib::Node(0,0,0.1));

    std::array<MeshLib::Node*, 3> nodes_array = {{nodes[0], nodes[1], nodes[2]}};
    std::vector<MeshLib::Element*> elements;
    MeshLib::Element* elem(new MeshLib::Tri(nodes_array));
    elements.push_back(elem);
    MeshLib::Mesh mesh ("testmesh", nodes, elements);

    MeshLib::MeshRevision rev(mesh);
    MeshLib::Mesh* result = rev.simplifyMesh("new_mesh", 0.2);

    ASSERT_EQ(MeshLib::MeshElemType::LINE, result->getElement(0)->getGeomType());
    ASSERT_EQ(1, result->getElement(0)->getContent());
    ASSERT_EQ(2u, result->getNumberOfNodes());
    delete result;

    result = rev.simplifyMesh("new_mesh", 0.0999);
    ASSERT_EQ(MeshLib::MeshElemType::TRIANGLE, result->getElement(0)->getGeomType());
    ASSERT_EQ(0.05, result->getElement(0)->getContent());
    ASSERT_EQ(3u, result->getNumberOfNodes());
    delete result;
}

TEST(MeshEditing, NonPlanarQuad)
{
    std::vector<MeshLib::Node*> nodes;
    nodes.push_back(new MeshLib::Node(0,0,0));
    nodes.push_back(new MeshLib::Node(0,1,0));
    nodes.push_back(new MeshLib::Node(1,1,0.1));
    nodes.push_back(new MeshLib::Node(1,0,0));

    std::array<MeshLib::Node*, 4> nodes_array = {{nodes[0], nodes[1], nodes[2], nodes[3]}};
    std::vector<MeshLib::Element*> elements;
    MeshLib::Element* elem (new MeshLib::Quad(nodes_array));
    elements.push_back(elem);
    MeshLib::Mesh mesh("testmesh", nodes, elements);
    MeshLib::MeshRevision rev(mesh);
    MeshLib::Mesh* result = rev.simplifyMesh("new_mesh", 0.2);

    ASSERT_EQ(2u, result->getNumberOfElements());
    ASSERT_EQ(MeshLib::MeshElemType::TRIANGLE, result->getElement(1)->getGeomType());

    delete result;
}

TEST(MeshEditing, Quad2Line)
{
    std::vector<MeshLib::Node*> nodes;
    nodes.push_back(new MeshLib::Node(1,0,0));
    nodes.push_back(new MeshLib::Node(0,1,0));
    nodes.push_back(new MeshLib::Node(0,1,0.1));
    nodes.push_back(new MeshLib::Node(1,0,0.1));

    std::array<MeshLib::Node*, 4> nodes_array = {{nodes[0], nodes[1], nodes[2], nodes[3]}};
    std::vector<MeshLib::Element*> elements;
    MeshLib::Element* elem (new MeshLib::Quad(nodes_array));
    elements.push_back(elem);
    MeshLib::Mesh mesh("testmesh", nodes, elements);

    MeshLib::MeshRevision rev(mesh);
    MeshLib::Mesh* result = rev.simplifyMesh("new_mesh", 0.2);

    ASSERT_EQ(MeshLib::MeshElemType::LINE, result->getElement(0)->getGeomType());
    ASSERT_NEAR(1.414213562373095, result->getElement(0)->getContent(), std::numeric_limits<double>::epsilon());
    ASSERT_EQ(2u, result->getNumberOfNodes());

    delete result;
}

TEST(MeshEditing, Quad2Tri)
{
    std::vector<MeshLib::Node*> nodes;
    nodes.push_back(new MeshLib::Node(1,0,0));
    nodes.push_back(new MeshLib::Node(0,1,0));
    nodes.push_back(new MeshLib::Node(0,1,0.1));
    nodes.push_back(new MeshLib::Node(1,1,0.1));

    std::array<MeshLib::Node*, 4> nodes_array = {{nodes[0], nodes[1], nodes[2], nodes[3]}};
    std::vector<MeshLib::Element*> elements;
    MeshLib::Element* elem (new MeshLib::Quad(nodes_array));
    elements.push_back(elem);
    MeshLib::Mesh mesh("testmesh", nodes, elements);

    MeshLib::MeshRevision rev(mesh);
    MeshLib::Mesh* result = rev.simplifyMesh("new_mesh", 0.2);

    ASSERT_EQ(MeshLib::MeshElemType::TRIANGLE, result->getElement(0)->getGeomType());
    ASSERT_NEAR(0.5049752469181039, result->getElement(0)->getContent(), std::numeric_limits<double>::epsilon());
    ASSERT_EQ(3u, result->getNumberOfNodes());

    delete result;
}

TEST(MeshEditing, NonPlanarHex)
{
    std::vector<MeshLib::Node*> nodes;
    nodes.push_back(new MeshLib::Node(0,0,-0.5));
    nodes.push_back(new MeshLib::Node(1,0,0));
    nodes.push_back(new MeshLib::Node(1,1,0));
    nodes.push_back(new MeshLib::Node(0,1,0));
    nodes.push_back(new MeshLib::Node(0,0,1));
    nodes.push_back(new MeshLib::Node(1,0,1));
    nodes.push_back(new MeshLib::Node(1,1,1));
    nodes.push_back(new MeshLib::Node(0,1,1));

    std::array<MeshLib::Node*, 8> nodes_array = {{nodes[0], nodes[1], nodes[2], nodes[3], nodes[4], nodes[5], nodes[6], nodes[7]}};
    std::vector<MeshLib::Element*> elements;
    MeshLib::Element* elem (new MeshLib::Hex(nodes_array));
    elements.push_back(elem);
    MeshLib::Mesh mesh("testmesh", nodes, elements);

    MeshLib::MeshRevision rev(mesh);
    MeshLib::Mesh* result = rev.simplifyMesh("new_mesh", 0.2);

    ASSERT_EQ(6u, result->getNumberOfElements());
    ASSERT_EQ(MeshLib::MeshElemType::TETRAHEDRON, result->getElement(4)->getGeomType());
    ASSERT_NEAR(0.25, result->getElement(0)->getContent(), std::numeric_limits<double>::epsilon());
    ASSERT_NEAR(0.1666666666666667, result->getElement(5)->getContent(), std::numeric_limits<double>::epsilon());

    delete result;
}

TEST(MeshEditing, Hex2PyramidPrism)
{
    std::vector<MeshLib::Node*> nodes;
    nodes.push_back(new MeshLib::Node(0,0,0));
    nodes.push_back(new MeshLib::Node(1,0,0));
    nodes.push_back(new MeshLib::Node(1,1,0));
    nodes.push_back(new MeshLib::Node(0,1,0));
    nodes.push_back(new MeshLib::Node(0,0,1));
    nodes.push_back(new MeshLib::Node(1,0,.1));
    nodes.push_back(new MeshLib::Node(1,1,1));
    nodes.push_back(new MeshLib::Node(0,1,1));

    std::array<MeshLib::Node*, 8> nodes_array = {{nodes[0], nodes[1], nodes[2], nodes[3], nodes[4], nodes[5], nodes[6], nodes[7]}};
    std::vector<MeshLib::Element*> elements;
    MeshLib::Element* elem (new MeshLib::Hex(nodes_array));
    elements.push_back(elem);
    MeshLib::Mesh mesh("testmesh", nodes, elements);

    MeshLib::MeshRevision rev(mesh);
    MeshLib::Mesh* result = rev.simplifyMesh("new_mesh", 0.2);

    ASSERT_EQ(2u, result->getNumberOfElements());
    ASSERT_EQ(MeshLib::MeshElemType::PYRAMID, result->getElement(0)->getGeomType());
    ASSERT_EQ(MeshLib::MeshElemType::PRISM, result->getElement(1)->getGeomType());
    ASSERT_NEAR(0.3333333333333333, result->getElement(0)->getContent(), std::numeric_limits<double>::epsilon());
    ASSERT_NEAR(0.5, result->getElement(1)->getContent(), std::numeric_limits<double>::epsilon());

    delete result;
}

TEST(MeshEditing, Hex2FourTets)
{
    std::vector<MeshLib::Node*> nodes;
    nodes.push_back(new MeshLib::Node(0,0,0));
    nodes.push_back(new MeshLib::Node(1,1,0));
    nodes.push_back(new MeshLib::Node(1,1,0));
    nodes.push_back(new MeshLib::Node(0,1,0));
    nodes.push_back(new MeshLib::Node(1,0,1));
    nodes.push_back(new MeshLib::Node(1,0,1));
    nodes.push_back(new MeshLib::Node(1,1,1));
    nodes.push_back(new MeshLib::Node(0,1,1));

    std::array<MeshLib::Node*, 8> nodes_array = {{nodes[0], nodes[1], nodes[2], nodes[3], nodes[4], nodes[5], nodes[6], nodes[7]}};
    std::vector<MeshLib::Element*> elements;
    MeshLib::Element* elem (new MeshLib::Hex(nodes_array));
    elements.push_back(elem);
    MeshLib::Mesh mesh("testmesh", nodes, elements);

    MeshLib::MeshRevision rev(mesh);
    MeshLib::Mesh* result = rev.simplifyMesh("new_mesh", 0.2);

    ASSERT_EQ(4u, result->getNumberOfElements());
    ASSERT_EQ(MeshLib::MeshElemType::TETRAHEDRON, result->getElement(1)->getGeomType());
    ASSERT_NEAR(0.1666666666666667, result->getElement(0)->getContent(), std::numeric_limits<double>::epsilon());
    ASSERT_NEAR(0.1666666666666667, result->getElement(1)->getContent(), std::numeric_limits<double>::epsilon());
    ASSERT_NEAR(0.1666666666666667, result->getElement(2)->getContent(), std::numeric_limits<double>::epsilon());
    ASSERT_NEAR(0.1666666666666667, result->getElement(3)->getContent(), std::numeric_limits<double>::epsilon());

    delete result;
}

TEST(MeshEditing, Hex2TwoTets)
{
    std::vector<MeshLib::Node*> nodes;
    nodes.push_back(new MeshLib::Node(0,0,0));
    nodes.push_back(new MeshLib::Node(1,0,0));
    nodes.push_back(new MeshLib::Node(1,1,0));
    nodes.push_back(new MeshLib::Node(0,1,0));
    nodes.push_back(new MeshLib::Node(1,1,1));
    nodes.push_back(new MeshLib::Node(1,1,1));
    nodes.push_back(new MeshLib::Node(1,1,1));
    nodes.push_back(new MeshLib::Node(1,1,1));

    std::array<MeshLib::Node*, 8> nodes_array = {{nodes[0], nodes[1], nodes[2], nodes[3], nodes[4], nodes[5], nodes[6], nodes[7]}};
    std::vector<MeshLib::Element*> elements;
    MeshLib::Element* elem (new MeshLib::Hex(nodes_array));
    elements.push_back(elem);
    MeshLib::Mesh mesh("testmesh", nodes, elements);

    MeshLib::MeshRevision rev(mesh);
    MeshLib::Mesh* result = rev.simplifyMesh("new_mesh", 0.2);

    ASSERT_EQ(2u, result->getNumberOfElements());
    ASSERT_EQ(MeshLib::MeshElemType::TETRAHEDRON, result->getElement(1)->getGeomType());
    ASSERT_NEAR(0.1666666666666667, result->getElement(0)->getContent(), std::numeric_limits<double>::epsilon());
    ASSERT_NEAR(0.1666666666666667, result->getElement(1)->getContent(), std::numeric_limits<double>::epsilon());

    delete result;
}

TEST(MeshEditing, NonPlanarPyramid)
{
    std::vector<MeshLib::Node*> nodes;
    nodes.push_back(new MeshLib::Node(0,0,0));
    nodes.push_back(new MeshLib::Node(1,0,-.5));
    nodes.push_back(new MeshLib::Node(1,1,0));
    nodes.push_back(new MeshLib::Node(0,1,0));
    nodes.push_back(new MeshLib::Node(1,0,1));

    std::array<MeshLib::Node*, 5> nodes_array = {{nodes[0], nodes[1], nodes[2], nodes[3], nodes[4]}};
    std::vector<MeshLib::Element*> elements;
    MeshLib::Element* elem (new MeshLib::Pyramid(nodes_array));
    elements.push_back(elem);
    MeshLib::Mesh mesh("testmesh", nodes, elements);

    MeshLib::MeshRevision rev(mesh);
    MeshLib::Mesh* result = rev.simplifyMesh("new_mesh", 0.2);

    ASSERT_EQ(2u, result->getNumberOfElements());
    ASSERT_EQ(MeshLib::MeshElemType::TETRAHEDRON, result->getElement(1)->getGeomType());
    ASSERT_NEAR(0.25, result->getElement(0)->getContent(),  std::numeric_limits<double>::epsilon());
    ASSERT_NEAR(0.1666666666666667, result->getElement(1)->getContent(), std::numeric_limits<double>::epsilon());

    delete result;
}

TEST(MeshEditing, Pyramid2Tet)
{
    std::vector<MeshLib::Node*> nodes;
    nodes.push_back(new MeshLib::Node(1,0,0));
    nodes.push_back(new MeshLib::Node(0,1,0));
    nodes.push_back(new MeshLib::Node(1,1,0));
    nodes.push_back(new MeshLib::Node(1,0,0));
    nodes.push_back(new MeshLib::Node(1,0,1));

    std::array<MeshLib::Node*, 5> nodes_array = {{nodes[0], nodes[1], nodes[2], nodes[3], nodes[4]}};
    std::vector<MeshLib::Element*> elements;
    MeshLib::Element* elem (new MeshLib::Pyramid(nodes_array));
    elements.push_back(elem);
    MeshLib::Mesh mesh("testmesh", nodes, elements);

    MeshLib::MeshRevision rev(mesh);
    MeshLib::Mesh* result = rev.simplifyMesh("new_mesh", 0.2);

    ASSERT_EQ(1u, result->getNumberOfElements());
    ASSERT_EQ(MeshLib::MeshElemType::TETRAHEDRON, result->getElement(0)->getGeomType());
    ASSERT_NEAR(0.16666666666666666, result->getElement(0)->getContent(), std::numeric_limits<double>::epsilon());

    delete result;
}

TEST(MeshEditing, Pyramid2Quad)
{
    std::vector<MeshLib::Node*> nodes;
    nodes.push_back(new MeshLib::Node(1,0,0));
    nodes.push_back(new MeshLib::Node(1,1,0));
    nodes.push_back(new MeshLib::Node(0,1,0));
    nodes.push_back(new MeshLib::Node(0,0,0));
    nodes.push_back(new MeshLib::Node(1,0,0.1));

    std::array<MeshLib::Node*, 5> nodes_array = {{nodes[0], nodes[1], nodes[2], nodes[3], nodes[4]}};
    std::vector<MeshLib::Element*> elements;
    MeshLib::Element* elem (new MeshLib::Pyramid(nodes_array));
    elements.push_back(elem);
    MeshLib::Mesh mesh("testmesh", nodes, elements);

    MeshLib::MeshRevision rev(mesh);
    MeshLib::Mesh* result = rev.simplifyMesh("new_mesh", 0.2);

    ASSERT_EQ(1u, result->getNumberOfElements());
    ASSERT_EQ(MeshLib::MeshElemType::QUAD, result->getElement(0)->getGeomType());
    ASSERT_NEAR(1, result->getElement(0)->getContent(), std::numeric_limits<double>::epsilon());

    delete result;
}

TEST(MeshEditing, Pyramid2Tri)
{
    std::vector<MeshLib::Node*> nodes;
    nodes.push_back(new MeshLib::Node(1,0,0));
    nodes.push_back(new MeshLib::Node(1,0.1,0));
    nodes.push_back(new MeshLib::Node(0,1,0));
    nodes.push_back(new MeshLib::Node(0,0,0));
    nodes.push_back(new MeshLib::Node(1,0,0.1));

    std::array<MeshLib::Node*, 5> nodes_array = {{nodes[0], nodes[1], nodes[2], nodes[3], nodes[4]}};
    std::vector<MeshLib::Element*> elements;
    MeshLib::Element* elem (new MeshLib::Pyramid(nodes_array));
    elements.push_back(elem);
    MeshLib::Mesh mesh("testmesh", nodes, elements);

    MeshLib::MeshRevision rev(mesh);
    MeshLib::Mesh* result = rev.simplifyMesh("new_mesh", 0.2);

    ASSERT_EQ(3u, result->getNumberOfNodes());
    ASSERT_EQ(1u, result->getNumberOfElements());
    ASSERT_EQ(MeshLib::MeshElemType::TRIANGLE, result->getElement(0)->getGeomType());
    ASSERT_NEAR(0.5, result->getElement(0)->getContent(), std::numeric_limits<double>::epsilon());

    delete result;
}

TEST(MeshEditing, NonPlanarPrism)
{
    std::vector<MeshLib::Node*> nodes;
    nodes.push_back(new MeshLib::Node(1,0,0));
    nodes.push_back(new MeshLib::Node(1,1,0));
    nodes.push_back(new MeshLib::Node(0,0,0));
    nodes.push_back(new MeshLib::Node(1,0,1));
    nodes.push_back(new MeshLib::Node(1,1,1));
    nodes.push_back(new MeshLib::Node(0,-0.5,2));

    std::array<MeshLib::Node*, 6> nodes_array = {{nodes[0], nodes[1], nodes[2], nodes[3], nodes[4], nodes[5]}};
    std::vector<MeshLib::Element*> elements;
    MeshLib::Element* elem (new MeshLib::Prism(nodes_array));
    elements.push_back(elem);
    MeshLib::Mesh mesh("testmesh", nodes, elements);

    MeshLib::MeshRevision rev(mesh);
    MeshLib::Mesh* result = rev.simplifyMesh("new_mesh", 0.2);

    ASSERT_EQ(3u, result->getNumberOfElements());
    ASSERT_EQ(MeshLib::MeshElemType::TETRAHEDRON, result->getElement(2)->getGeomType());
    ASSERT_NEAR(0.1666666666666667, result->getElement(0)->getContent(), std::numeric_limits<double>::epsilon());

    delete result;
}

TEST(MeshEditing, Prism2TwoTets)
{
    std::vector<MeshLib::Node*> nodes;
    nodes.push_back(new MeshLib::Node(1,0,0));
    nodes.push_back(new MeshLib::Node(1,1,0));
    nodes.push_back(new MeshLib::Node(0,0,0));
    nodes.push_back(new MeshLib::Node(1,0.9,1));
    nodes.push_back(new MeshLib::Node(1,1,1));
    nodes.push_back(new MeshLib::Node(0,0,1));

    std::array<MeshLib::Node*, 6> nodes_array = {{nodes[0], nodes[1], nodes[2], nodes[3], nodes[4], nodes[5]}};
    std::vector<MeshLib::Element*> elements;
    MeshLib::Element* elem (new MeshLib::Prism(nodes_array));
    elements.push_back(elem);
    MeshLib::Mesh mesh("testmesh", nodes, elements);

    MeshLib::MeshRevision rev(mesh);
    MeshLib::Mesh* result = rev.simplifyMesh("new_mesh", 0.2);

    ASSERT_EQ(5u, result->getNumberOfNodes());
    ASSERT_EQ(2u, result->getNumberOfElements());
    ASSERT_EQ(MeshLib::MeshElemType::TETRAHEDRON, result->getElement(1)->getGeomType());
    ASSERT_NEAR(0.1666666666666667, result->getElement(0)->getContent(), std::numeric_limits<double>::epsilon());
    ASSERT_NEAR(0.15, result->getElement(1)->getContent(), std::numeric_limits<double>::epsilon());

    delete result;
}

TEST(MeshEditing, Prism2Quad)
{
    std::vector<MeshLib::Node*> nodes;
    nodes.push_back(new MeshLib::Node(1,0.9,0));
    nodes.push_back(new MeshLib::Node(1,1,0));
    nodes.push_back(new MeshLib::Node(0,0,0));
    nodes.push_back(new MeshLib::Node(1,0.9,1));
    nodes.push_back(new MeshLib::Node(1,1,1));
    nodes.push_back(new MeshLib::Node(0,0,1));

    std::array<MeshLib::Node*, 6> nodes_array = {{nodes[0], nodes[1], nodes[2], nodes[3], nodes[4], nodes[5]}};
    std::vector<MeshLib::Element*> elements;
    MeshLib::Element* elem (new MeshLib::Prism(nodes_array));
    elements.push_back(elem);
    MeshLib::Mesh mesh("testmesh", nodes, elements);

    MeshLib::MeshRevision rev(mesh);
    MeshLib::Mesh* result = rev.simplifyMesh("new_mesh", 0.2);

    ASSERT_EQ(1u, result->getNumberOfElements());
    ASSERT_EQ(MeshLib::MeshElemType::QUAD, result->getElement(0)->getGeomType());
    ASSERT_NEAR(1.345362404707371, result->getElement(0)->getContent(), std::numeric_limits<double>::epsilon());

    delete result;
}

TEST(MeshEditing, Prism2Tet)
{
    std::vector<MeshLib::Node*> nodes;
    nodes.push_back(new MeshLib::Node(1,0,0));
    nodes.push_back(new MeshLib::Node(1,1,0));
    nodes.push_back(new MeshLib::Node(0.9,0,0));
    nodes.push_back(new MeshLib::Node(1,0.9,1));
    nodes.push_back(new MeshLib::Node(1,1,1));
    nodes.push_back(new MeshLib::Node(0,0,1));

    std::array<MeshLib::Node*, 6> nodes_array = {{nodes[0], nodes[1], nodes[2], nodes[3], nodes[4], nodes[5]}};
    std::vector<MeshLib::Element*> elements;
    MeshLib::Element* elem (new MeshLib::Prism(nodes_array));
    elements.push_back(elem);
    MeshLib::Mesh mesh("testmesh", nodes, elements);

    MeshLib::MeshRevision rev(mesh);
    MeshLib::Mesh* result = rev.simplifyMesh("new_mesh", 0.2);

    ASSERT_EQ(1u, result->getNumberOfElements());
    ASSERT_EQ(MeshLib::MeshElemType::TETRAHEDRON, result->getElement(0)->getGeomType());
    ASSERT_NEAR(0.1666666666666667, result->getElement(0)->getContent(), std::numeric_limits<double>::epsilon());

    delete result;
}

TEST(MeshEditing, Prism2Tri)
{
    std::vector<MeshLib::Node*> nodes;
    nodes.push_back(new MeshLib::Node(1,0,0));
    nodes.push_back(new MeshLib::Node(1,1,0));
    nodes.push_back(new MeshLib::Node(1,1,0));
    nodes.push_back(new MeshLib::Node(1,1,1));
    nodes.push_back(new MeshLib::Node(1,1,1));
    nodes.push_back(new MeshLib::Node(0.9,0.9,1));

    std::array<MeshLib::Node*, 6> nodes_array = {{nodes[0], nodes[1], nodes[2], nodes[3], nodes[4], nodes[5]}};
    std::vector<MeshLib::Element*> elements;
    MeshLib::Element* elem (new MeshLib::Prism(nodes_array));
    elements.push_back(elem);
    MeshLib::Mesh mesh("testmesh", nodes, elements);

    MeshLib::MeshRevision rev(mesh);
    MeshLib::Mesh* result = rev.simplifyMesh("new_mesh", 0.2);

    ASSERT_EQ(1u, result->getNumberOfElements());
    ASSERT_EQ(MeshLib::MeshElemType::TRIANGLE, result->getElement(0)->getGeomType());
    ASSERT_NEAR(0.5, result->getElement(0)->getContent(), std::numeric_limits<double>::epsilon());

    delete result;
}
