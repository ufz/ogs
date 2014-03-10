/**
 * @file TestMeshRevision.cpp
 * @author Karsten Rink
 * @date 2013-03-04
 * @brief Tests for MeshRevision class
 *
 * @copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include "gtest/gtest.h"

#include "Mesh.h"
#include "Node.h"
#include "MeshEditing/MeshRevision.h"
#include "Elements/Element.h"
#include "Elements/Tri.h"
#include "Elements/Quad.h"
#include "Elements/Tet.h"
#include "Elements/Hex.h"
#include "Elements/Pyramid.h"
#include "Elements/Prism.h"


TEST(MeshEditing, Tri)
{
	std::vector<MeshLib::Node*> nodes;
	nodes.push_back(new MeshLib::Node(1,0,0));
	nodes.push_back(new MeshLib::Node(0,0,0));
	nodes.push_back(new MeshLib::Node(0,0,0.1));

	MeshLib::Node* tri_nodes[3] = { nodes[0], nodes[1], nodes[2] };
	std::vector<MeshLib::Element*> elements;
	MeshLib::Element* elem (new MeshLib::Tri(tri_nodes));
	elements.push_back(elem);
	MeshLib::Mesh mesh ("testmesh", nodes, elements);

	MeshLib::MeshRevision rev(mesh);
	MeshLib::Mesh* result = rev.simplifyMesh("new_mesh", 0.2);

	ASSERT_EQ(result->getElement(0)->getGeomType(), MeshElemType::LINE);
	ASSERT_EQ(result->getElement(0)->getContent(), 1);
	ASSERT_EQ(result->getNNodes(), 2);

	delete result;
}

TEST(MeshEditing, NonPlanarQuad)
{
	std::vector<MeshLib::Node*> nodes;
	nodes.push_back(new MeshLib::Node(0,0,0));
	nodes.push_back(new MeshLib::Node(0,1,0));
	nodes.push_back(new MeshLib::Node(1,1,0.1));
	nodes.push_back(new MeshLib::Node(1,0,0));

	MeshLib::Node* nodes_array[4] = { nodes[0], nodes[1], nodes[2], nodes[3] };
	std::vector<MeshLib::Element*> elements;
	MeshLib::Element* elem (new MeshLib::Quad(nodes_array));
	elements.push_back(elem);
	MeshLib::Mesh mesh("testmesh", nodes, elements);
	MeshLib::MeshRevision rev(mesh);
	MeshLib::Mesh* result = rev.simplifyMesh("new_mesh", 0.2);

	ASSERT_EQ(result->getNElements(), 2);
	ASSERT_EQ(result->getElement(1)->getGeomType(), MeshElemType::TRIANGLE);

	delete result;
}

TEST(MeshEditing, Quad2Line)
{
	std::vector<MeshLib::Node*> nodes;
	nodes.push_back(new MeshLib::Node(1,0,0));
	nodes.push_back(new MeshLib::Node(0,1,0));
	nodes.push_back(new MeshLib::Node(0,1,0.1));
	nodes.push_back(new MeshLib::Node(1,0,0.1));

	MeshLib::Node* nodes_array[4] = { nodes[0], nodes[1], nodes[2], nodes[3] };
	std::vector<MeshLib::Element*> elements;
	MeshLib::Element* elem (new MeshLib::Quad(nodes_array));
	elements.push_back(elem);
	MeshLib::Mesh mesh("testmesh", nodes, elements);

	MeshLib::MeshRevision rev(mesh);
	MeshLib::Mesh* result = rev.simplifyMesh("new_mesh", 0.2);
	
	ASSERT_EQ(result->getElement(0)->getGeomType(), MeshElemType::LINE);
	ASSERT_NEAR(result->getElement(0)->getContent(), 1.414213562373095, std::numeric_limits<double>::epsilon());
	ASSERT_EQ(result->getNNodes(), 2);

	delete result;
}

TEST(MeshEditing, Quad2Tri)
{
	std::vector<MeshLib::Node*> nodes;
	nodes.push_back(new MeshLib::Node(1,0,0));
	nodes.push_back(new MeshLib::Node(0,1,0));
	nodes.push_back(new MeshLib::Node(0,1,0.1));
	nodes.push_back(new MeshLib::Node(1,1,0.1));

	MeshLib::Node* nodes_array[4] = { nodes[0], nodes[1], nodes[2], nodes[3] };
	std::vector<MeshLib::Element*> elements;
	MeshLib::Element* elem (new MeshLib::Quad(nodes_array));
	elements.push_back(elem);
	MeshLib::Mesh mesh("testmesh", nodes, elements);

	MeshLib::MeshRevision rev(mesh);
	MeshLib::Mesh* result = rev.simplifyMesh("new_mesh", 0.2);

	ASSERT_EQ(result->getElement(0)->getGeomType(), MeshElemType::TRIANGLE);
	ASSERT_NEAR(result->getElement(0)->getContent(), 0.5049752469181039, std::numeric_limits<double>::epsilon());
	ASSERT_EQ(result->getNNodes(), 3);

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

	MeshLib::Node* nodes_array[8] = { nodes[0], nodes[1], nodes[2], nodes[3], nodes[4], nodes[5], nodes[6], nodes[7] };
	std::vector<MeshLib::Element*> elements;
	MeshLib::Element* elem (new MeshLib::Hex(nodes_array));
	elements.push_back(elem);
	MeshLib::Mesh mesh("testmesh", nodes, elements);

	MeshLib::MeshRevision rev(mesh);
	MeshLib::Mesh* result = rev.simplifyMesh("new_mesh", 0.2);

	ASSERT_EQ(result->getNElements(), 6);
	ASSERT_EQ(result->getElement(4)->getGeomType(), MeshElemType::TETRAHEDRON);
	ASSERT_NEAR(result->getElement(0)->getContent(), 0.25, std::numeric_limits<double>::epsilon());
	ASSERT_NEAR(result->getElement(5)->getContent(), 0.1666666666666667, std::numeric_limits<double>::epsilon());

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

	MeshLib::Node* nodes_array[8] = { nodes[0], nodes[1], nodes[2], nodes[3], nodes[4], nodes[5], nodes[6], nodes[7] };
	std::vector<MeshLib::Element*> elements;
	MeshLib::Element* elem (new MeshLib::Hex(nodes_array));
	elements.push_back(elem);
	MeshLib::Mesh mesh("testmesh", nodes, elements);

	MeshLib::MeshRevision rev(mesh);
	MeshLib::Mesh* result = rev.simplifyMesh("new_mesh", 0.2);

	ASSERT_EQ(result->getNElements(), 2);
	ASSERT_EQ(result->getElement(0)->getGeomType(), MeshElemType::PYRAMID);
	ASSERT_EQ(result->getElement(1)->getGeomType(), MeshElemType::PRISM);
	ASSERT_NEAR(result->getElement(0)->getContent(), 0.3333333333333333, std::numeric_limits<double>::epsilon());
	ASSERT_NEAR(result->getElement(1)->getContent(), 0.5, std::numeric_limits<double>::epsilon());
	
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

	MeshLib::Node* nodes_array[8] = { nodes[0], nodes[1], nodes[2], nodes[3], nodes[4], nodes[5], nodes[6], nodes[7] };
	std::vector<MeshLib::Element*> elements;
	MeshLib::Element* elem (new MeshLib::Hex(nodes_array));
	elements.push_back(elem);
	MeshLib::Mesh mesh("testmesh", nodes, elements);

	MeshLib::MeshRevision rev(mesh);
	MeshLib::Mesh* result = rev.simplifyMesh("new_mesh", 0.2);

	ASSERT_EQ(result->getNElements(), 4);
	ASSERT_EQ(result->getElement(1)->getGeomType(), MeshElemType::TETRAHEDRON);
	ASSERT_NEAR(result->getElement(0)->getContent(), 0.1666666666666667, std::numeric_limits<double>::epsilon());
	ASSERT_NEAR(result->getElement(1)->getContent(), 0.1666666666666667, std::numeric_limits<double>::epsilon());
	ASSERT_NEAR(result->getElement(2)->getContent(), 0.1666666666666667, std::numeric_limits<double>::epsilon());
	ASSERT_NEAR(result->getElement(3)->getContent(), 0.1666666666666667, std::numeric_limits<double>::epsilon());
	
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

	MeshLib::Node* nodes_array[8] = { nodes[0], nodes[1], nodes[2], nodes[3], nodes[4], nodes[5], nodes[6], nodes[7] };
	std::vector<MeshLib::Element*> elements;
	MeshLib::Element* elem (new MeshLib::Hex(nodes_array));
	elements.push_back(elem);
	MeshLib::Mesh mesh("testmesh", nodes, elements);

	MeshLib::MeshRevision rev(mesh);
	MeshLib::Mesh* result = rev.simplifyMesh("new_mesh", 0.2);

	ASSERT_EQ(result->getNElements(), 2);
	ASSERT_EQ(result->getElement(1)->getGeomType(), MeshElemType::TETRAHEDRON);
	ASSERT_NEAR(result->getElement(0)->getContent(), 0.1666666666666667, std::numeric_limits<double>::epsilon());
	ASSERT_NEAR(result->getElement(1)->getContent(), 0.1666666666666667, std::numeric_limits<double>::epsilon());
	
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

	MeshLib::Node* nodes_array[5] = { nodes[0], nodes[1], nodes[2], nodes[3], nodes[4] };
	std::vector<MeshLib::Element*> elements;
	MeshLib::Element* elem (new MeshLib::Pyramid(nodes_array));
	elements.push_back(elem);
	MeshLib::Mesh mesh("testmesh", nodes, elements);

	MeshLib::MeshRevision rev(mesh);
	MeshLib::Mesh* result = rev.simplifyMesh("new_mesh", 0.2);

	ASSERT_EQ(result->getNElements(), 2);
	ASSERT_EQ(result->getElement(1)->getGeomType(), MeshElemType::TETRAHEDRON);
	ASSERT_NEAR(result->getElement(0)->getContent(),  0.25, std::numeric_limits<double>::epsilon());
	ASSERT_NEAR(result->getElement(1)->getContent(), 0.1666666666666667, std::numeric_limits<double>::epsilon());
	
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

	MeshLib::Node* nodes_array[5] = { nodes[0], nodes[1], nodes[2], nodes[3], nodes[4] };
	std::vector<MeshLib::Element*> elements;
	MeshLib::Element* elem (new MeshLib::Pyramid(nodes_array));
	elements.push_back(elem);
	MeshLib::Mesh mesh("testmesh", nodes, elements);

	MeshLib::MeshRevision rev(mesh);
	MeshLib::Mesh* result = rev.simplifyMesh("new_mesh", 0.2);

	ASSERT_EQ(result->getNElements(), 1);
	ASSERT_EQ(result->getElement(0)->getGeomType(), MeshElemType::TETRAHEDRON);
	ASSERT_NEAR(result->getElement(0)->getContent(), 0.16666666666666666, std::numeric_limits<double>::epsilon());
	
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

	MeshLib::Node* nodes_array[5] = { nodes[0], nodes[1], nodes[2], nodes[3], nodes[4] };
	std::vector<MeshLib::Element*> elements;
	MeshLib::Element* elem (new MeshLib::Pyramid(nodes_array));
	elements.push_back(elem);
	MeshLib::Mesh mesh("testmesh", nodes, elements);

	MeshLib::MeshRevision rev(mesh);
	MeshLib::Mesh* result = rev.simplifyMesh("new_mesh", 0.2);

	ASSERT_EQ(result->getNElements(), 1);
	ASSERT_EQ(result->getElement(0)->getGeomType(), MeshElemType::QUAD);
	ASSERT_NEAR(result->getElement(0)->getContent(), 1, std::numeric_limits<double>::epsilon());
	
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

	MeshLib::Node* nodes_array[5] = { nodes[0], nodes[1], nodes[2], nodes[3], nodes[4] };
	std::vector<MeshLib::Element*> elements;
	MeshLib::Element* elem (new MeshLib::Pyramid(nodes_array));
	elements.push_back(elem);
	MeshLib::Mesh mesh("testmesh", nodes, elements);

	MeshLib::MeshRevision rev(mesh);
	MeshLib::Mesh* result = rev.simplifyMesh("new_mesh", 0.2);

	ASSERT_EQ(result->getNNodes(), 3);
	ASSERT_EQ(result->getNElements(), 1);
	ASSERT_EQ(result->getElement(0)->getGeomType(), MeshElemType::TRIANGLE);
	ASSERT_NEAR(result->getElement(0)->getContent(), 0.5, std::numeric_limits<double>::epsilon());
	
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

	MeshLib::Node* nodes_array[6] = { nodes[0], nodes[1], nodes[2], nodes[3], nodes[4], nodes[5] };
	std::vector<MeshLib::Element*> elements;
	MeshLib::Element* elem (new MeshLib::Prism(nodes_array));
	elements.push_back(elem);
	MeshLib::Mesh mesh("testmesh", nodes, elements);

	MeshLib::MeshRevision rev(mesh);
	MeshLib::Mesh* result = rev.simplifyMesh("new_mesh", 0.2);

	ASSERT_EQ(result->getNElements(), 3);
	ASSERT_EQ(result->getElement(2)->getGeomType(), MeshElemType::TETRAHEDRON);
	ASSERT_NEAR(result->getElement(0)->getContent(), 0.1666666666666667, std::numeric_limits<double>::epsilon());
	
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

	MeshLib::Node* nodes_array[6] = { nodes[0], nodes[1], nodes[2], nodes[3], nodes[4], nodes[5] };
	std::vector<MeshLib::Element*> elements;
	MeshLib::Element* elem (new MeshLib::Prism(nodes_array));
	elements.push_back(elem);
	MeshLib::Mesh mesh("testmesh", nodes, elements);

	MeshLib::MeshRevision rev(mesh);
	MeshLib::Mesh* result = rev.simplifyMesh("new_mesh", 0.2);

	ASSERT_EQ(result->getNNodes(), 5);
	ASSERT_EQ(result->getNElements(), 2);
	ASSERT_EQ(result->getElement(1)->getGeomType(), MeshElemType::TETRAHEDRON);
	ASSERT_NEAR(result->getElement(0)->getContent(), 0.1666666666666667, std::numeric_limits<double>::epsilon());
	ASSERT_NEAR(result->getElement(1)->getContent(), 0.15, std::numeric_limits<double>::epsilon());
	
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

	MeshLib::Node* nodes_array[6] = { nodes[0], nodes[1], nodes[2], nodes[3], nodes[4], nodes[5] };
	std::vector<MeshLib::Element*> elements;
	MeshLib::Element* elem (new MeshLib::Prism(nodes_array));
	elements.push_back(elem);
	MeshLib::Mesh mesh("testmesh", nodes, elements);

	MeshLib::MeshRevision rev(mesh);
	MeshLib::Mesh* result = rev.simplifyMesh("new_mesh", 0.2);

	ASSERT_EQ(result->getNElements(), 1);
	ASSERT_EQ(result->getElement(0)->getGeomType(), MeshElemType::QUAD);
	ASSERT_NEAR(result->getElement(0)->getContent(), 1.345362404707371, std::numeric_limits<double>::epsilon());	
	
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

	MeshLib::Node* nodes_array[6] = { nodes[0], nodes[1], nodes[2], nodes[3], nodes[4], nodes[5] };
	std::vector<MeshLib::Element*> elements;
	MeshLib::Element* elem (new MeshLib::Prism(nodes_array));
	elements.push_back(elem);
	MeshLib::Mesh mesh("testmesh", nodes, elements);

	MeshLib::MeshRevision rev(mesh);
	MeshLib::Mesh* result = rev.simplifyMesh("new_mesh", 0.2);

	ASSERT_EQ(result->getNElements(), 1);
	ASSERT_EQ(result->getElement(0)->getGeomType(), MeshElemType::TETRAHEDRON);
	ASSERT_NEAR(result->getElement(0)->getContent(), 0.1666666666666667, std::numeric_limits<double>::epsilon());
	
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

	MeshLib::Node* nodes_array[6] = { nodes[0], nodes[1], nodes[2], nodes[3], nodes[4], nodes[5] };
	std::vector<MeshLib::Element*> elements;
	MeshLib::Element* elem (new MeshLib::Prism(nodes_array));
	elements.push_back(elem);
	MeshLib::Mesh mesh("testmesh", nodes, elements);

	MeshLib::MeshRevision rev(mesh);
	MeshLib::Mesh* result = rev.simplifyMesh("new_mesh", 0.2);

	ASSERT_EQ(result->getNElements(), 1);
	ASSERT_EQ(result->getElement(0)->getGeomType(), MeshElemType::TRIANGLE);
	ASSERT_NEAR(result->getElement(0)->getContent(), 0.5, std::numeric_limits<double>::epsilon());

	delete result;
}
