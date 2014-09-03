/**
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include "gtest/gtest.h"

#include "Mesh.h"
#include "MeshSubsets.h"

using namespace MeshLib;

TEST(MeshLibMeshSubsets, UniqueMeshIds)
{
	// Create first mesh
	Mesh const m0("first", std::vector<Node*>(), std::vector<Element*>());
	Mesh const m1("second", std::vector<Node*>(), std::vector<Element*>());

	std::vector<Node*> const empty_node_ptr_vector;

	MeshSubset const ms0(m0, empty_node_ptr_vector);
	MeshSubset const ms1(m1, empty_node_ptr_vector);
	MeshSubset const ms1a(m1, empty_node_ptr_vector);

	MeshSubset const* const mesh_subsets[3] = {&ms0, &ms1, &ms1a};

	EXPECT_NO_THROW(MeshSubsets(&mesh_subsets[0], &mesh_subsets[0] + 2));
	EXPECT_THROW(MeshSubsets(&mesh_subsets[1], &mesh_subsets[1] + 2), std::logic_error);
	EXPECT_THROW(MeshSubsets(&mesh_subsets[0], &mesh_subsets[0] + 3), std::logic_error);
}

