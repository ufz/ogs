/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <vtkCellType.h>
#include <vtkCellTypeSource.h>
#include <vtkNew.h>
#include <vtkTransform.h>
#include <vtkTransformFilter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>

#include "MeshLib/Elements/Element.h"
#include "MeshLib/Elements/TemplateElement.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshGenerators/VtkMeshConverter.h"
#include "Tests/Utils.h"

inline VTKCellType toVtk(MeshLib::CellType const cell_type)
{
    using namespace MeshLib;

    switch (cell_type)
    {
        case CellType::LINE2:
            return VTK_LINE;
        case CellType::TRI3:
            return VTK_TRIANGLE;
        case CellType::QUAD4:
            return VTK_QUAD;
        case CellType::TET4:
            return VTK_TETRA;
        case CellType::HEX8:
            return VTK_HEXAHEDRON;
        case CellType::PRISM6:
            return VTK_WEDGE;
        case CellType::PYRAMID5:
            return VTK_PYRAMID;
        default:
            OGS_FATAL("Unsupported cell type " + CellType2String(cell_type));
    }
}

template <typename ElementRule>
MeshLib::CellType getCellType(Type<MeshLib::TemplateElement<ElementRule>>)
{
    return ElementRule::cell_type;
}

inline void checkBounds(vtkUnstructuredGrid& grid, std::size_t const dim)
{
    double bounds[6];
    grid.GetBounds(bounds);

    for (std::size_t c = 0; c < dim; ++c)
    {
        auto const min = bounds[2 * c];
        auto const max = bounds[2 * c + 1];

        if (min != -0.5 || max != 0.5)
        {
            OGS_FATAL("Unexpected bounds for a unit cube. " +
                      std::to_string(min) + " != -0.5 or " +
                      std::to_string(max) + " != 0.5 for component " +
                      std::to_string(c) + '.');
        }
    }
}

inline void checkCellTypes(MeshLib::Mesh const& mesh,
                           MeshLib::CellType const& cell_type)
{
    for (std::size_t i = 0; i < mesh.getNumberOfElements(); ++i)
    {
        auto const& e = *mesh.getElement(i);
        if (e.getCellType() != cell_type)
        {
            OGS_FATAL("Unexpected element type: " +
                      MeshLib::CellType2String(e.getCellType()) +
                      " != " + MeshLib::CellType2String(cell_type));
        }
    }
}

inline void checkVolume(MeshLib::Mesh const& mesh)
{
    double volume = 0;
    for (std::size_t i = 0; i < mesh.getNumberOfElements(); ++i)
    {
        auto const& e = *mesh.getElement(i);
        volume += e.getContent();
    }

    auto const tol = std::numeric_limits<double>::epsilon();
    auto const diff = std::abs(1.0 - volume);
    if (diff > tol)
    {
        OGS_FATAL(
            "The volume of a unit cube must be 1. Instead, the volume is " +
            std::to_string(volume) + ". The difference is " +
            std::to_string(diff) + " > " + std::to_string(tol) + '.');
    }
}

// Creates a 1/2/3 dimensional unit cube centered at (0, 0, 0) subdivided into
// one or more elements of the given MeshElementType.
template <typename MeshElementType>
std::unique_ptr<MeshLib::Mesh> createUnitCube()
{
    auto const cell_type = getCellType(Type<MeshElementType>{});
    auto const dim = MeshElementType::dimension;

    // See
    // https://kitware.github.io/vtk-examples/site/Cxx/GeometricObjects/CellTypeSource/
    vtkNew<vtkCellTypeSource> source;
    source->SetCellType(toVtk(cell_type));

    vtkNew<vtkTransform> transform;
    transform->Translate(-0.5, dim >= 2 ? -0.5 : 0.0, dim >= 3 ? -0.5 : 0.0);

    vtkNew<vtkTransformFilter> transform_filter;
    transform_filter->SetInputConnection(source->GetOutputPort());
    transform_filter->SetTransform(transform);
    transform_filter->Update();

    auto* vtk_grid =
        dynamic_cast<vtkUnstructuredGrid*>(transform_filter->GetOutput());

    checkBounds(*vtk_grid, dim);

    std::unique_ptr<MeshLib::Mesh> ogs_grid{
        MeshLib::VtkMeshConverter::convertUnstructuredGrid(vtk_grid,
                                                           "unit_cube")};

    checkCellTypes(*ogs_grid, cell_type);
    checkVolume(*ogs_grid);

    return ogs_grid;
}
