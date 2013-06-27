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

#include <vector>
#include <gtest/gtest.h>

#include "MeshLib/MeshGenerators/MeshGenerator.h"

#include "VecMatOnMeshLib/VecMeshItems/VectorComposition.h"
#include "VecMatOnMeshLib/VecMeshItems/ComponentDistribution.h"

TEST(VecMatOnMeshLib, DataArrangementByComponentType)
{
	// This test checks the following case:
	// - Each node has two scalar values (i.e. two components per node)
	// - A vector is created such that it contains all nodal values
	// - Entries in the vector are arranged in the order of a component type and then node ID
	// - For example, x=[(node 0, comp 0) (node 1, comp 0) ... (node n, comp0), (node 0, comp1) ... ]

    //assume the following mesh having line elements
    std::unique_ptr<MeshLib::Mesh> msh(MeshLib::MeshGenerator::generateLineMesh(1.0, 9));

    //data component 0 and 1 are assigned to all nodes in the mesh
    const std::size_t comp0_id = 0;
    const std::size_t comp1_id = 1;
    const VecMatOnMeshLib::MeshItems mesh_items(msh.get(), msh->getNodes());
    VecMatOnMeshLib::ComponentDistribution comp0(&mesh_items);
    VecMatOnMeshLib::ComponentDistribution comp1(&mesh_items);

    //set up data arrangement
    std::vector<VecMatOnMeshLib::ComponentDistribution*> vec_comp_dis;
    vec_comp_dis.push_back(&comp0);
    vec_comp_dis.push_back(&comp1);
    VecMatOnMeshLib::VectorComposition da(vec_comp_dis, VecMatOnMeshLib::OrderingType::BY_COMPONENT_TYPE);
    //std::cout << "# database \n";
    //da.print();
    std::cout << std::endl;

    //check
    ASSERT_EQ(20u, da.size());
    ASSERT_EQ(2u, da.getNComponents());
    ASSERT_EQ(1u, da.getNMeshes());
    ASSERT_EQ(0u, da.getDataID(VecMatOnMeshLib::MeshItem(msh->getID(), VecMatOnMeshLib::MeshItemType::Node, 0), comp0_id));
    ASSERT_EQ(1u, da.getDataID(VecMatOnMeshLib::MeshItem(msh->getID(), VecMatOnMeshLib::MeshItemType::Node, 1), comp0_id));
    ASSERT_EQ(10u, da.getDataID(VecMatOnMeshLib::MeshItem(msh->getID(), VecMatOnMeshLib::MeshItemType::Node, 0), comp1_id));
    ASSERT_EQ(11u, da.getDataID(VecMatOnMeshLib::MeshItem(msh->getID(), VecMatOnMeshLib::MeshItemType::Node, 1), comp1_id));
    auto vecCompIDs = da.getComponentIDs(VecMatOnMeshLib::MeshItem(msh->getID(), VecMatOnMeshLib::MeshItemType::Node, 0));
    ASSERT_EQ(2u, vecCompIDs.size());
    ASSERT_EQ(0u, vecCompIDs[0]);
    ASSERT_EQ(1u, vecCompIDs[1]);

    //check out of range
    ASSERT_EQ(std::numeric_limits<std::size_t>::max(), da.getDataID(VecMatOnMeshLib::MeshItem(msh->getID(), VecMatOnMeshLib::MeshItemType::Node, 10), comp0_id));
    ASSERT_EQ(std::numeric_limits<std::size_t>::max(), da.getDataID(VecMatOnMeshLib::MeshItem(msh->getID()+1, VecMatOnMeshLib::MeshItemType::Node, 0), comp0_id));
    ASSERT_EQ(std::numeric_limits<std::size_t>::max(), da.getDataID(VecMatOnMeshLib::MeshItem(msh->getID(), VecMatOnMeshLib::MeshItemType::Cell, 0), comp0_id));
    ASSERT_EQ(std::numeric_limits<std::size_t>::max(), da.getDataID(VecMatOnMeshLib::MeshItem(msh->getID(), VecMatOnMeshLib::MeshItemType::Node, 0), 10));
}

TEST(VecMatOnMeshLib, DataArrangementByMeshItem)
{
	// This test checks the following case:
	// - Each node has two scalar values (i.e. two components per node)
	// - A vector is created such that it contains all nodal values
	// - Entries in the vector are arranged in the order of node ID and then a component type
	// - For example, x=[(node 0, comp 0) (node 0, comp 1) ... (node n, comp0), (node n, comp1) ]

    //assume the following mesh
    std::unique_ptr<MeshLib::Mesh> msh(MeshLib::MeshGenerator::generateLineMesh(1.0, 9));

    //data component 0 and 1 are assigned to all nodes in the mesh
    const std::size_t comp0_id = 0;
    const std::size_t comp1_id = 1;
    const VecMatOnMeshLib::MeshItems mesh_items(msh.get(), msh->getNodes());
    VecMatOnMeshLib::ComponentDistribution comp0(&mesh_items);
    VecMatOnMeshLib::ComponentDistribution comp1(&mesh_items);

    //set up data arrangement
    std::vector<VecMatOnMeshLib::ComponentDistribution*> vec_comp_dis;
    vec_comp_dis.push_back(&comp0);
    vec_comp_dis.push_back(&comp1);
    VecMatOnMeshLib::VectorComposition da(vec_comp_dis, VecMatOnMeshLib::OrderingType::BY_MESH_ITEM_ID);
    //std::cout << "# database \n";
    //da.print();
    std::cout << std::endl;

    //check
    ASSERT_EQ(20u, da.size());
    ASSERT_EQ(2u, da.getNComponents());
    ASSERT_EQ(1u, da.getNMeshes());
    ASSERT_EQ(0u, da.getDataID(VecMatOnMeshLib::MeshItem(msh->getID(), VecMatOnMeshLib::MeshItemType::Node, 0), comp0_id));
    ASSERT_EQ(1u, da.getDataID(VecMatOnMeshLib::MeshItem(msh->getID(), VecMatOnMeshLib::MeshItemType::Node, 0), comp1_id));
    ASSERT_EQ(2u, da.getDataID(VecMatOnMeshLib::MeshItem(msh->getID(), VecMatOnMeshLib::MeshItemType::Node, 1), comp0_id));
    ASSERT_EQ(3u, da.getDataID(VecMatOnMeshLib::MeshItem(msh->getID(), VecMatOnMeshLib::MeshItemType::Node, 1), comp1_id));
}

