/**
 * @copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include "gtest/gtest.h"

#include <memory>
#include <cmath>

#include "Mesh.h"
#include "Node.h"
#include "MeshGenerators/MeshGenerator.h"
#include "Elements/Element.h"
#include "Elements/Hex.h"

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
}

