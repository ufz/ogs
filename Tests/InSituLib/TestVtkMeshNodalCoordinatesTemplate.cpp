/**
 * \file
 * \author Lars Bilke
 * \date   2014-02-26
 * \brief  Unit tests for In-Situ mesh nodes
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */
#include <vtkNew.h>

#include "gtest/gtest.h"

#include "Mesh.h"
#include "MeshGenerators/MeshGenerator.h"

#include "VtkMeshNodalCoordinatesTemplate.h"

TEST(InSituLibNodalCoordinates, Init)
{
	MeshLib::Mesh* mesh = MeshLib::MeshGenerator::generateRegularQuadMesh(10.0, 99);

	vtkNew<VtkMeshNodalCoordinatesTemplate<double> > nodeCoords;
	nodeCoords->SetNodes(mesh->getNodes());
	nodeCoords->PrintSelf(std::cout, vtkIndent());

	delete mesh;
}
