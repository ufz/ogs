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

using namespace MeshLib;

TEST(MeshLib, ElementConstatns)
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

