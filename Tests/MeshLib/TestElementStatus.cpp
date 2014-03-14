/**
 * @file TestElementStatus.cpp
 * @author Karsten Rink
 * @date 2013-03-14
 * @brief Tests for ElementStatus class
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
#include "Elements/Element.h"
#include "ElementStatus.h"
#include "MeshGenerators/MeshGenerator.h"


TEST(MeshLib, ElementStatus)
{
const unsigned width (100);
	const unsigned elements_per_side (20);
	MeshLib::Mesh* mesh (MeshLib::MeshGenerator::generateRegularQuadMesh(width, elements_per_side));
	MeshLib::ElementStatus status(mesh);
	const std::vector<MeshLib::Element*> elements (mesh->getElements());

	for (unsigned i=0; i<elements_per_side; ++i)
	{
		for (unsigned j=0; j<elements_per_side; ++j)
			elements[i*elements_per_side + j]->setValue(i);
	}

	// all elements active
	ASSERT_EQ (status.getNActiveElements(), elements.size());
	// all nodes active
	ASSERT_EQ (status.getNActiveNodes(), mesh->getNNodes());

	// set material 1 to false
	status.setMaterialStatus(1, false);
	ASSERT_EQ (status.getNActiveElements(), elements.size()-elements_per_side);

	// set material 1 to false (again)
	status.setMaterialStatus(1, false);
	ASSERT_EQ (status.getNActiveElements(), elements.size()-elements_per_side);

	// set material 0 to false
	status.setMaterialStatus(0, false);
	ASSERT_EQ (status.getNActiveElements(), elements.size()-(2*elements_per_side));

	// active elements
	std::vector<unsigned> active_elements (status.getActiveElements());
	ASSERT_EQ (active_elements.size(), status.getNActiveElements());

	// active nodes
	std::vector<unsigned> active_nodes (status.getActiveNodes());
	ASSERT_EQ (active_nodes.size(), status.getNActiveNodes());

	// set element 1 to false (yet again)
	status.setElementStatus(1, false);
	status.isActive(1);
	ASSERT_EQ (status.getNActiveElements(), elements.size()-(2*elements_per_side));
	ASSERT_EQ (status.getNActiveNodes(), mesh->getNNodes()-(2*(elements_per_side+1)));

	// set element 1 to true
	status.setElementStatus(1, true);
	ASSERT_EQ (status.getNActiveElements(), elements.size()-(2*elements_per_side)+1);
	ASSERT_EQ (status.getNActiveNodes(), mesh->getNNodes()-(2*(elements_per_side+1))+4);
	ASSERT_EQ(status.isActive(1), true);

	std::vector<unsigned> active_elements_at_node (status.getActiveElementsAtNode(2));
	ASSERT_EQ(active_elements_at_node.size(), 1);
	active_elements_at_node = status.getActiveElementsAtNode(22);
	ASSERT_EQ(active_elements_at_node.size(), 1);
	active_elements_at_node = status.getActiveElementsAtNode(102);
	ASSERT_EQ(active_elements_at_node.size(), 4);

	status.setAll(true);
	ASSERT_EQ(status.getNActiveElements(), elements.size());
	ASSERT_EQ(status.getNActiveNodes(), mesh->getNNodes());

	status.setAll(false);
	ASSERT_EQ(0, status.getNActiveElements());
	ASSERT_EQ(0, status.getNActiveNodes());
}
