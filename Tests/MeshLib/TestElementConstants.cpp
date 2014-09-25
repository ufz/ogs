/**
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include "gtest/gtest.h"

#include "Elements/Quad.h"

using namespace MeshLib;

TEST(MeshLib, Quad4ElementConstants)
{
	ASSERT_EQ(2u, Quad::dimension);
	ASSERT_EQ(4u, Quad::n_all_nodes);
	ASSERT_EQ(4u, Quad::n_base_nodes);

	std::array<Node*, 4> nodes;
	nodes[0] = new Node(0.0, 0.0, 0.0);
	nodes[1] = new Node(0.0, 1.0, 0.0);
	nodes[2] = new Node(1.0, 1.0, 0.0);
	nodes[3] = new Node(1.0, 0.0, 0.0);
	Quad quad(nodes);

	ASSERT_EQ(2u, quad.getDimension());
	ASSERT_EQ(4u, quad.getNNodes());
	ASSERT_EQ(4u, quad.getNLinearNodes());
}

TEST(MeshLib, Quad8ElementConstants)
{
	ASSERT_EQ(2u, Quad8::dimension);
	ASSERT_EQ(8u, Quad8::n_all_nodes);
	ASSERT_EQ(4u, Quad8::n_base_nodes);

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
	ASSERT_EQ(8u, quad8.getNNodes());
	ASSERT_EQ(4u, quad8.getNLinearNodes());
}

TEST(MeshLib, Quad9ElementConstants)
{
	ASSERT_EQ(2u, Quad9::dimension);
	ASSERT_EQ(9u, Quad9::n_all_nodes);
	ASSERT_EQ(4u, Quad9::n_base_nodes);

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
	ASSERT_EQ(9u, quad9.getNNodes());
	ASSERT_EQ(4u, quad9.getNLinearNodes());
}

