/**
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include <ctime>
#include "gtest/gtest.h"

#include "Elements/Quad.h"
#include "Elements/Tri.h"
#include "MshEnums.h"

using namespace MeshLib;

TEST(MeshLib, Quad4ElementConstants)
{
	ASSERT_EQ(2u, Quad::_dimension);
	ASSERT_EQ(4u, Quad::_n_all_nodes);
	ASSERT_EQ(4u, Quad::_n_base_nodes);

	std::array<Node*, 4> nodes;
	nodes[0] = new Node(0.0, 0.0, 0.0);
	nodes[1] = new Node(0.0, 1.0, 0.0);
	nodes[2] = new Node(1.0, 1.0, 0.0);
	nodes[3] = new Node(1.0, 0.0, 0.0);
	Quad quad(nodes);

	ASSERT_EQ(2u, quad.getDimension());
	ASSERT_EQ(4u, quad.getNNodes(true));
	ASSERT_EQ(4u, quad.getNNodes(false));
}

TEST(MeshLib, Quad8ElementConstants)
{
	typedef TemplateQuad<8, CellType::QUAD8> Quad8;
	ASSERT_EQ(2u, Quad8::_dimension);
	ASSERT_EQ(8u, Quad8::_n_all_nodes);
	ASSERT_EQ(4u, Quad8::_n_base_nodes);

	std::array<Node*, 8> nodes;
	nodes[0] = new Node(0.0, 0.0, 0.0);
	nodes[1] = new Node(0.0, 1.0, 0.0);
	nodes[2] = new Node(1.0, 1.0, 0.0);
	nodes[3] = new Node(1.0, 0.0, 0.0);

	nodes[4] = new Node(0.5, 0.0, 0.0);
	nodes[5] = new Node(1.0, 0.5, 0.0);
	nodes[6] = new Node(0.5, 1.0, 0.0);
	nodes[7] = new Node(0.5, 0.0, 0.0);

	Quad8 quad8(nodes);

	ASSERT_EQ(2u, quad8.getDimension());
	ASSERT_EQ(8u, quad8.getNNodes(true));
	ASSERT_EQ(4u, quad8.getNNodes(false));
}

TEST(MeshLib, Quad9ElementConstants)
{
	typedef TemplateQuad<9, CellType::QUAD9> Quad9;
	ASSERT_EQ(2u, Quad9::_dimension);
	ASSERT_EQ(9u, Quad9::_n_all_nodes);
	ASSERT_EQ(4u, Quad9::_n_base_nodes);

	std::array<Node*, 9> nodes;
	nodes[0] = new Node(0.0, 0.0, 0.0);
	nodes[1] = new Node(0.0, 1.0, 0.0);
	nodes[2] = new Node(1.0, 1.0, 0.0);
	nodes[3] = new Node(1.0, 0.0, 0.0);

	nodes[4] = new Node(0.5, 0.0, 0.0);
	nodes[5] = new Node(1.0, 0.5, 0.0);
	nodes[6] = new Node(0.5, 1.0, 0.0);
	nodes[7] = new Node(0.5, 0.0, 0.0);

	nodes[8] = new Node(0.5, 0.5, 0.0);
	Quad9 quad9(nodes);

	ASSERT_EQ(2u, quad9.getDimension());
	ASSERT_EQ(9u, quad9.getNNodes(true));
	ASSERT_EQ(4u, quad9.getNNodes(false));
}

