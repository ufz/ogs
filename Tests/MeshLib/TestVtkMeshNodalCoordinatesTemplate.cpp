/**
 * \file
 * \author Lars Bilke
 * \date   2014-02-26
 * \brief  Unit tests for In-Situ mesh nodes
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */
#include <vtkNew.h>
#include <vtkUnstructuredGrid.h>

#include "gtest/gtest.h"

#include "MeshLib/Mesh.h"
#include "MeshLib/MeshGenerators/MeshGenerator.h"
#include "MeshLib/Vtk/VtkMeshNodalCoordinatesTemplate.h"

TEST(MeshLibdalCoordinates, Init)
{
    const std::size_t subdivisions = 99;
    const double length = 10.0;
    const double dx = length / subdivisions;

    MeshLib::Mesh* mesh = MeshLib::MeshGenerator::generateRegularQuadMesh(length, subdivisions);

    vtkNew<MeshLib::VtkMeshNodalCoordinatesTemplate<double> > nodeCoords;
    nodeCoords->SetNodes(mesh->getNodes());
    //nodeCoords->PrintSelf(std::cout, vtkIndent());

    ASSERT_EQ(nodeCoords->GetNumberOfComponents(), 3);
    const std::size_t numTuples = (subdivisions+1)*(subdivisions+1);
    ASSERT_EQ(nodeCoords->GetNumberOfTuples(), numTuples);

    // First point
    ASSERT_EQ(nodeCoords->GetValueReference(0), 0.0);
    ASSERT_EQ(nodeCoords->GetValueReference(1), 0.0);
    ASSERT_EQ(nodeCoords->GetValueReference(2), 0.0);

    // First second
    ASSERT_NEAR(nodeCoords->GetValueReference(3), dx, std::numeric_limits<double>::epsilon());
    ASSERT_NEAR(nodeCoords->GetValueReference(4), 0.0, std::numeric_limits<double>::epsilon());
    ASSERT_NEAR(nodeCoords->GetValueReference(5), 0.0, std::numeric_limits<double>::epsilon());

    double* coords = nodeCoords->GetTuple(1);
    ASSERT_NEAR(coords[0], dx, std::numeric_limits<double>::epsilon());
    ASSERT_NEAR(coords[1], 0.0, std::numeric_limits<double>::epsilon());
    ASSERT_NEAR(coords[2], 0.0, std::numeric_limits<double>::epsilon());

    delete mesh;
}
