/**
 * @file TestElementStatus.cpp
 * @author Karsten Rink
 * @date 2013-03-14
 * @brief Tests for ElementStatus class
 *
 * @copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include <memory>

#include "gtest/gtest.h"

#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/ElementStatus.h"
#include "MeshLib/MeshGenerators/MeshGenerator.h"


TEST(MeshLib, ElementStatus)
{
	const unsigned width (100);
	const unsigned elements_per_side (20);
	auto const mesh = std::unique_ptr<MeshLib::Mesh>{
		MeshLib::MeshGenerator::generateRegularQuadMesh(width, elements_per_side)};

	boost::optional<MeshLib::PropertyVector<int> &> material_id_properties(
		mesh->getProperties().createNewPropertyVector<int>("MaterialIDs",
		                                                   MeshLib::MeshItemType::Cell)
	);
	ASSERT_FALSE(!material_id_properties);
	(*material_id_properties).resize(mesh->getNElements());

	const std::vector<MeshLib::Element*> elements (mesh->getElements());

	for (unsigned i=0; i<elements_per_side; ++i)
	{
		for (unsigned j=0; j<elements_per_side; ++j)
			(*material_id_properties)[elements[i*elements_per_side + j]->getID()] = i;
	}

	{
		// all elements and nodes active
		MeshLib::ElementStatus status(mesh.get());
		ASSERT_EQ (elements.size(), status.getNActiveElements());
		ASSERT_EQ (mesh->getNNodes(), status.getNActiveNodes());
	}

	{
		// set material 1 to false
		std::vector<unsigned> inactiveMat{1};
		MeshLib::ElementStatus status(mesh.get(), inactiveMat);
		ASSERT_EQ (elements.size()-elements_per_side, status.getNActiveElements());
	}

	{
		// set material 0 and 1 to false
		std::vector<unsigned> inactiveMat{0, 1};
		MeshLib::ElementStatus status(mesh.get(), inactiveMat);
		ASSERT_EQ (elements.size()-(2*elements_per_side), status.getNActiveElements());

		// active elements
		auto &active_elements (status.getActiveElements());
		ASSERT_EQ (active_elements.size(), status.getNActiveElements());

		// active nodes
		auto& active_nodes (status.getActiveNodes());
		ASSERT_EQ (active_nodes.size(), status.getNActiveNodes());
	}
}
