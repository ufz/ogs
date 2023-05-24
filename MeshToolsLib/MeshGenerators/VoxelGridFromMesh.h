/**
 * \file
 * \copyright
 * Copyright (c) 2012-2023, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <vtkAbstractArray.h>
#include <vtkCellData.h>
#include <vtkCellLocator.h>
#include <vtkDataArray.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridReader.h>

#include "InfoLib/GitInfo.h"
#include "MeshLib/Mesh.h"
#include "MeshToolsLib/MeshEditing/RemoveMeshComponents.h"
#include "MeshToolsLib/MeshGenerators/MeshGenerator.h"
#include "MeshLib/MeshSearch/ElementSearch.h"
#include "MeshLib/MeshSearch/MeshElementGrid.h"

static constexpr std::string const cell_id_name = "CellIds";

namespace MeshToolsLib::MeshGenerator::VoxelGridFromMesh
{

    std::array<std::size_t, 3> getDimensions(
        MathLib::Point3d const& min,
        MathLib::Point3d const& max,
        std::array<double, 3> const& cellsize);

    std::vector<int> assignCellIds(
        vtkSmartPointer<vtkUnstructuredGrid> const& mesh,
        MathLib::Point3d const& min,
        std::array<std::size_t, 3> const& dims,
        std::array<double, 3> const& cellsize);

    bool removeUnusedGridCells(vtkSmartPointer<vtkUnstructuredGrid> const& mesh,
                               std::unique_ptr<MeshLib::Mesh>& grid);

    void mapMeshArraysOntoGrid(vtkSmartPointer<vtkUnstructuredGrid> const& mesh,
                               std::unique_ptr<MeshLib::Mesh>& grid);
};