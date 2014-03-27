/**
 * @file TestDuplicate.cpp
 * @author Karsten Rink
 * @date 2013-03-25
 * @brief Tests for Duplicate functions
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
#include "MeshEditing/DuplicateMeshComponents.h"
#include "MeshGenerators/MeshGenerator.h"
#include "MeshQuality/MeshValidation.h"
#include "MathTools.h"

TEST(MeshLib, Duplicate)
{
	MeshLib::Mesh* mesh (MeshLib::MeshGenerator::generateRegularQuadMesh(10, 5, 1));

	std::vector<MeshLib::Node*> new_nodes (MeshLib::copyNodeVector(mesh->getNodes()));
	std::vector<MeshLib::Element*> new_elements (MeshLib::copyElementVector(mesh->getElements(), new_nodes));

	MeshLib::Mesh new_mesh ("new", new_nodes, new_elements);

	ASSERT_EQ (mesh->getNElements(), new_mesh.getNElements());
	ASSERT_EQ (mesh->getNNodes(), new_mesh.getNNodes());

	std::vector<std::size_t> del_idx(1,1);
	MeshLib::removeMeshNodes(*mesh, del_idx);

	ASSERT_EQ (mesh->getNElements(), new_mesh.getNElements()-2);
	ASSERT_EQ (mesh->getNNodes(), new_mesh.getNNodes()-2);

	ASSERT_DOUBLE_EQ (4.0, MathLib::sqrDist(*mesh->getNode(0), *new_mesh.getNode(0)));
	ASSERT_DOUBLE_EQ (0.0, MathLib::sqrDist(*mesh->getNode(0), *new_mesh.getNode(2)));

	ASSERT_DOUBLE_EQ (4.0, MathLib::sqrDist(*mesh->getElement(0)->getNode(0), *new_mesh.getElement(0)->getNode(0)));
	ASSERT_DOUBLE_EQ (0.0, MathLib::sqrDist(*mesh->getElement(0)->getNode(0), *new_mesh.getElement(2)->getNode(0)));
}
