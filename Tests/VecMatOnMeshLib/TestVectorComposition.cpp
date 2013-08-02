/**
 * \file
 * \author Norihiro Watanabe
 * \date   2013-04-16
 * \brief  Implementation tests.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <gtest/gtest.h>
#include <vector>

#include "MeshLib/MeshGenerators/MeshGenerator.h"

#include "VecMatOnMeshLib/VecMeshItems/MeshComponentMap.h"
#include "VecMatOnMeshLib/VecMeshItems/MeshSubsets.h"

class VecMatOnMeshLibTest : public ::testing::Test
{
	public:
	VecMatOnMeshLibTest()
	    : mesh(nullptr), mesh_items(nullptr)
	{
		mesh = MeshLib::MeshGenerator::generateLineMesh(1.0, 9);
		mesh_items = new VecMatOnMeshLib::MeshSubset(*mesh, mesh->getNodes());

		//set up data arrangement
		vec_comp_dis.emplace_back(new VecMatOnMeshLib::MeshSubsets(mesh_items));
		vec_comp_dis.emplace_back(new VecMatOnMeshLib::MeshSubsets(mesh_items));

	}
	~VecMatOnMeshLibTest()
	{
		std::remove_if(vec_comp_dis.begin(), vec_comp_dis.end(),
		    [](VecMatOnMeshLib::MeshSubsets* p) { delete p; return true; });
		delete mesh_items;
		delete mesh;
	}

	MeshLib::Mesh const* mesh;
	VecMatOnMeshLib::MeshSubset const* mesh_items;

	//data component 0 and 1 are assigned to all nodes in the mesh
	std::size_t const comp0_id = 0;
	std::size_t const comp1_id = 1;
	std::vector<VecMatOnMeshLib::MeshSubsets*> vec_comp_dis;
};

TEST_F(VecMatOnMeshLibTest, DataArrangementByComponentType)
{
	// This test checks the following case:
	// - Each node has two scalar values (i.e. two components per node)
	// - A vector is created such that it contains all nodal values
	// - Entries in the vector are arranged in the order of a component type and then node ID
	// - For example, x=[(node 0, comp 0) (node 1, comp 0) ... (node n, comp0), (node 0, comp1) ... ]

	//set up data arrangement
	VecMatOnMeshLib::MeshComponentMap da(
	    vec_comp_dis, VecMatOnMeshLib::ComponentOrder::BY_COMPONENT);
	//std::cout << "# database \n";
	//da.print();
	//std::cout << std::endl;

	using VecMatOnMeshLib::MeshItemType;
	typedef VecMatOnMeshLib::Location Loc;
	//check
	auto checkNodeAndComponent =
	    [&da, this](std::size_t const n, std::size_t const c)
	    {
	        return da.getDataID(Loc(mesh->getID(), MeshItemType::Node, n), c);
	    };

	ASSERT_EQ(20u, da.size());
	ASSERT_EQ(0u , checkNodeAndComponent(0, comp0_id));
	ASSERT_EQ(1u , checkNodeAndComponent(1, comp0_id));
	ASSERT_EQ(10u, checkNodeAndComponent(0, comp1_id));
	ASSERT_EQ(11u, checkNodeAndComponent(1, comp1_id));

	auto vecCompIDs = da.getComponentIDs(
	    Loc(mesh->getID(), MeshItemType::Node, 0));
	ASSERT_EQ(2u, vecCompIDs.size());
	ASSERT_EQ(0u, vecCompIDs[0]);
	ASSERT_EQ(1u, vecCompIDs[1]);

	//check out of range
	std::size_t const out_of_range = std::numeric_limits<std::size_t>::max();
	ASSERT_EQ(out_of_range, checkNodeAndComponent(10, comp0_id));
	ASSERT_EQ(out_of_range, da.getDataID(
	    Loc(mesh->getID() + 1, MeshItemType::Node, 0), comp0_id));
	ASSERT_EQ(out_of_range, da.getDataID(
	    Loc(mesh->getID(), MeshItemType::Cell, 0), comp0_id));
	ASSERT_EQ(out_of_range, da.getDataID(
	    Loc(mesh->getID(), MeshItemType::Node, 0), 10));
}

TEST_F(VecMatOnMeshLibTest, DataArrangementByLocationType)
{
	// This test checks the following case:
	// - Each node has two scalar values (i.e. two components per node)
	// - A vector is created such that it contains all nodal values
	// - Entries in the vector are arranged in the order of node ID and then a component type
	// - For example, x=[(node 0, comp 0) (node 0, comp 1) ... (node n, comp0), (node n, comp1) ]

	//set up data arrangement
	VecMatOnMeshLib::MeshComponentMap da(
	    vec_comp_dis, VecMatOnMeshLib::ComponentOrder::BY_LOCATION);
	//std::cout << "# database \n";
	//da.print();
	//std::cout << std::endl;

	using VecMatOnMeshLib::MeshItemType;
	typedef VecMatOnMeshLib::Location Loc;
	//check
	auto checkNodeAndComponent =
	    [&da, this](std::size_t const n, std::size_t const c)
	    {
	        return da.getDataID(Loc(mesh->getID(), MeshItemType::Node, n), c);
	    };

	ASSERT_EQ(20u, da.size());
	ASSERT_EQ(0u, checkNodeAndComponent(0, comp0_id));
	ASSERT_EQ(1u, checkNodeAndComponent(0, comp1_id));
	ASSERT_EQ(2u, checkNodeAndComponent(1, comp0_id));
	ASSERT_EQ(3u, checkNodeAndComponent(1, comp1_id));
	auto vecCompIDs = da.getComponentIDs(
	    Loc(mesh->getID(), MeshItemType::Node, 0));
	ASSERT_EQ(2u, vecCompIDs.size());
	ASSERT_EQ(0u, vecCompIDs[0]);
	ASSERT_EQ(1u, vecCompIDs[1]);

	//check out of range
	std::size_t const out_of_range = std::numeric_limits<std::size_t>::max();
	ASSERT_EQ(out_of_range, checkNodeAndComponent(10, comp0_id));
	ASSERT_EQ(out_of_range, da.getDataID(
	    Loc(mesh->getID() + 1, MeshItemType::Node, 0), comp0_id));
	ASSERT_EQ(out_of_range, da.getDataID(
	    Loc(mesh->getID(), MeshItemType::Cell, 0), comp0_id));
	ASSERT_EQ(out_of_range, da.getDataID(
	    Loc(mesh->getID(), MeshItemType::Node, 0), 10));
}
