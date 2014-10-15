/**
 * @copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include "gtest/gtest.h"

#include <memory>
#include <cmath>

#include "MathTools.h"
#include "Mesh.h"
#include "Node.h"
#include "MeshGenerators/MeshGenerator.h"
#include "Elements/Element.h"
#include "Elements/Hex.h"

#include "Tests/TestTools.h"

using namespace MeshLib;

TEST(MeshLib, MeshGeneratorRegularHex)
{
	const double L = 10.0;
	const std::size_t n_subdivisions = 9;
	const double dL = L / static_cast<double>(n_subdivisions);
	std::unique_ptr<Mesh> msh(MeshGenerator::generateRegularHexMesh(L, n_subdivisions));

	ASSERT_EQ(std::pow(n_subdivisions, 3), msh->getNElements());
	ASSERT_EQ(std::pow(n_subdivisions+1, 3), msh->getNNodes());

	// check nodes
	const Node& node0 = *msh->getNode(0);
	const Node& node1 = *msh->getNode(1);
	const Node& node_n = *msh->getNode(msh->getNNodes()-2);
	const Node& node_n1 = *msh->getNode(msh->getNNodes()-1);
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
	const Element& ele_n = *msh->getElement(msh->getNElements()-1);
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
	ASSERT_EQ(msh->getNNodes(), msh2->getNNodes());
	ASSERT_DOUBLE_EQ(0, MathLib::sqrDist(*(msh->getNode(msh->getNNodes()-1)), *(msh2->getNode(msh->getNNodes()-1))));

	unsigned n_x (10);
	unsigned n_y (5);
	unsigned n_z (2);
	double delta (1.2);
	std::unique_ptr<Mesh> hex_mesh (MeshGenerator::generateRegularHexMesh(n_x, n_y, n_z, delta));
	ASSERT_EQ(n_x * n_y * n_z, hex_mesh->getNElements());
	ASSERT_EQ((n_x+1) * (n_y+1) * (n_z+1), hex_mesh->getNNodes());
	const MeshLib::Node* node (hex_mesh->getNode(hex_mesh->getNNodes()-1));
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
	ASSERT_EQ(n_x * n_y, quad_mesh->getNElements());
	ASSERT_EQ((n_x+1) * (n_y+1), quad_mesh->getNNodes());
	const MeshLib::Node* node (quad_mesh->getNode(quad_mesh->getNNodes()-1));
	ASSERT_DOUBLE_EQ(n_x*delta, (*node)[0]);
	ASSERT_DOUBLE_EQ(n_y*delta, (*node)[1]);
	ASSERT_DOUBLE_EQ(0, (*node)[2]);

	const double L = 10.0;
	const std::size_t n_subdivisions = 9;
	std::unique_ptr<Mesh> quad_mesh2(MeshGenerator::generateRegularQuadMesh(L, n_subdivisions));
	ASSERT_EQ(n_subdivisions * n_subdivisions, quad_mesh2->getNElements());
	node = quad_mesh2->getNode(quad_mesh2->getNNodes()-1);
	ASSERT_DOUBLE_EQ(L, (*node)[0]);
	ASSERT_DOUBLE_EQ(L, (*node)[1]);
}

TEST(MeshLib, MeshGeneratorRegularQuad8)
{
	unsigned n_x (10);
	unsigned n_y (5);
	double delta (1.2);
	double tol(std::numeric_limits<float>::epsilon());
	std::unique_ptr<Mesh> quad_mesh (MeshGenerator::generateRegularQuad8Mesh(n_x, n_y,delta));
	ASSERT_EQ(n_x * n_y, quad_mesh->getNElements());
	ASSERT_EQ((n_x+1) * (n_y+1) + n_x * (n_y+1) + (n_x+1)*n_y, quad_mesh->getNNodes());
	const MeshLib::Element* ele0 = quad_mesh->getElement(0);
	ASSERT_EQ(CellType::QUAD8, ele0->getCellType());
	ASSERT_ARRAY_NEAR(Node(0,     0,     0), *ele0->getNode(0), 3, tol);
	ASSERT_ARRAY_NEAR(Node(delta, 0,     0), *ele0->getNode(1), 3, tol);
	ASSERT_ARRAY_NEAR(Node(delta, delta, 0), *ele0->getNode(2), 3, tol);
	ASSERT_ARRAY_NEAR(Node(0,     delta, 0), *ele0->getNode(3), 3, tol);
	ASSERT_ARRAY_NEAR(Node(delta*0.5, 0,         0), *ele0->getNode(4), 3, tol);
	ASSERT_ARRAY_NEAR(Node(delta,     delta*0.5, 0), *ele0->getNode(5), 3, tol);
	ASSERT_ARRAY_NEAR(Node(delta*0.5, delta,     0), *ele0->getNode(6), 3, tol);
	ASSERT_ARRAY_NEAR(Node(0,         delta*0.5, 0), *ele0->getNode(7), 3, tol);
	const MeshLib::Element* ele1 = quad_mesh->getElement(quad_mesh->getNElements()-1);
	ASSERT_EQ(CellType::QUAD8, ele1->getCellType());
	ASSERT_ARRAY_NEAR(Node(delta*(n_x-1), delta*(n_y-1), 0),         *ele1->getNode(0), 3, tol);
	ASSERT_ARRAY_NEAR(Node(delta*n_x,     delta*(n_y-1), 0),     *ele1->getNode(1), 3, tol);
	ASSERT_ARRAY_NEAR(Node(delta*n_x,     delta*n_y,     0), *ele1->getNode(2), 3, tol);
	ASSERT_ARRAY_NEAR(Node(delta*(n_x-1), delta*n_y,     0),     *ele1->getNode(3), 3, tol);
	ASSERT_ARRAY_NEAR(Node(delta*(n_x-0.5),  delta*(n_y-1),   0), *ele1->getNode(4), 3, tol);
	ASSERT_ARRAY_NEAR(Node(delta*n_x,        delta*(n_y-0.5), 0), *ele1->getNode(5), 3, tol);
	ASSERT_ARRAY_NEAR(Node(delta*(n_x-0.5),  delta*n_y,       0), *ele1->getNode(6), 3, tol);
	ASSERT_ARRAY_NEAR(Node(delta*(n_x-1),    delta*(n_y-0.5), 0), *ele1->getNode(7), 3, tol);
	const MeshLib::Node* nodel (quad_mesh->getNode((n_x+1)*(n_y+1)-1));
	ASSERT_DOUBLE_EQ(n_x*delta, (*nodel)[0]);
	ASSERT_DOUBLE_EQ(n_y*delta, (*nodel)[1]);
	ASSERT_DOUBLE_EQ(0, (*nodel)[2]);
	const MeshLib::Node* nodeq (quad_mesh->getNode(quad_mesh->getNNodes()-1));
	ASSERT_DOUBLE_EQ(n_x*delta-delta*0.5, (*nodeq)[0]);
	ASSERT_DOUBLE_EQ(n_y*delta, (*nodeq)[1]);
	ASSERT_DOUBLE_EQ(0, (*nodeq)[2]);

	const double L = 10.0;
	const std::size_t n_subdivisions = 9;
	std::unique_ptr<Mesh> quad_mesh2(MeshGenerator::generateRegularQuad8Mesh(L, n_subdivisions));
	ASSERT_EQ(n_subdivisions * n_subdivisions, quad_mesh2->getNElements());
	nodeq = quad_mesh2->getNode(quad_mesh2->getNNodes()-1);
	ASSERT_DOUBLE_EQ(L-L/n_subdivisions*0.5, (*nodeq)[0]);
	ASSERT_DOUBLE_EQ(L, (*nodeq)[1]);
	ASSERT_DOUBLE_EQ(0, (*nodeq)[2]);
}
