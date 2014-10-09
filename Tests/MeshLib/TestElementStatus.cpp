/**
 * @file TestElementStatus.cpp
 * @author Karsten Rink
 * @date 2013-03-14
 * @brief Tests for ElementStatus class
 *
 * @copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include "gtest/gtest.h"

#include "Mesh.h"
#include "MeshLib/Node.h"
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
	ASSERT_EQ (elements.size(), status.getNActiveElements());
	// all nodes active
	ASSERT_EQ (mesh->getNNodes(), status.getNActiveNodes());

	// set material 1 to false
	status.setMaterialStatus(1, false);
	ASSERT_EQ (elements.size()-elements_per_side, status.getNActiveElements());

	// set material 1 to false (again)
	status.setMaterialStatus(1, false);
	ASSERT_EQ (elements.size()-elements_per_side, status.getNActiveElements());

	// set material 0 to false
	status.setMaterialStatus(0, false);
	ASSERT_EQ (elements.size()-(2*elements_per_side), status.getNActiveElements());

	// active elements
	std::vector<std::size_t> active_elements (status.getActiveElements());
	ASSERT_EQ (active_elements.size(), status.getNActiveElements());

	// active nodes
	std::vector<std::size_t> active_nodes (status.getActiveNodes());
	ASSERT_EQ (active_nodes.size(), status.getNActiveNodes());

	// set element 1 to false (yet again)
	status.setElementStatus(1, false);
	status.getElementStatus(1);
	ASSERT_EQ (elements.size()-(2*elements_per_side), status.getNActiveElements());
	ASSERT_EQ (mesh->getNNodes()-(2*(elements_per_side+1)), status.getNActiveNodes());

	// set element 1 to true
	status.setElementStatus(1, true);
	ASSERT_EQ (elements.size()-(2*elements_per_side)+1, status.getNActiveElements());
	ASSERT_EQ (mesh->getNNodes()-(2*(elements_per_side+1))+4, status.getNActiveNodes());
	ASSERT_EQ(status.getElementStatus(1), true);

	std::vector<std::size_t> active_elements_at_node (status.getActiveElementsAtNode(2));
	ASSERT_EQ(1u, active_elements_at_node.size());
	active_elements_at_node = status.getActiveElementsAtNode(22);
	ASSERT_EQ(1u, active_elements_at_node.size());
	active_elements_at_node = status.getActiveElementsAtNode(44);
	ASSERT_EQ(2u, active_elements_at_node.size());
	active_elements_at_node = status.getActiveElementsAtNode(102);
	ASSERT_EQ(4u, active_elements_at_node.size());

	status.setAll(true);
	ASSERT_EQ(elements.size(), status.getNActiveElements());
	ASSERT_EQ(mesh->getNNodes(), status.getNActiveNodes());

	status.setAll(false);
	ASSERT_EQ(0u, status.getNActiveElements());
	ASSERT_EQ(0u, status.getNActiveNodes());
}
