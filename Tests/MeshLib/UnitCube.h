// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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
#include "MeshLib/IO/VtkIO/VtkMeshConverter.h"
#include "MeshLib/Mesh.h"
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
            OGS_FATAL("Unsupported cell type {}", CellType2String(cell_type));
    }
}

template <typename ElementRule>
MeshLib::CellType getCellType(
    Type<MeshLib::TemplateElement<ElementRule>> /*unused*/)
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
            OGS_FATAL(
                "Unexpected bounds for a unit cube. {} != -0.5 or {} != 0.5 "
                "for component {}.",
                min, max, c);
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
            OGS_FATAL("Unexpected element type: {} != {}",
                      MeshLib::CellType2String(e.getCellType()),
                      MeshLib::CellType2String(cell_type));
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
            "The volume of a unit cube must be 1. Instead, the volume is {}. "
            "The difference is {:g} > {:g}",
            volume, diff, tol);
    }
}

// Since VTK 9.7, vtkCellTypeSource generates wedges with the opposite node
// order of what VtkMeshConverter expects, i.e., inverted prisms (negative
// Jacobian determinant). Swap the triangle nodes of such prisms; prisms with
// the expected node order are left untouched, so this works with any VTK
// version.
inline void fixPrismNodeOrder(MeshLib::Mesh& mesh)
{
    for (auto* const e : mesh.getElements())
    {
        if (e->getCellType() != MeshLib::CellType::PRISM6)
        {
            continue;
        }

        auto const& n0 = *e->getNode(0);
        auto const& n1 = *e->getNode(1);
        auto const& n2 = *e->getNode(2);
        auto const& n3 = *e->getNode(3);

        // Triple product (n1 - n0) x (n2 - n0) . (n3 - n0).
        double const a[3] = {n1[0] - n0[0], n1[1] - n0[1], n1[2] - n0[2]};
        double const b[3] = {n2[0] - n0[0], n2[1] - n0[1], n2[2] - n0[2]};
        double const c[3] = {n3[0] - n0[0], n3[1] - n0[1], n3[2] - n0[2]};
        double const triple = (a[1] * b[2] - a[2] * b[1]) * c[0] +
                              (a[2] * b[0] - a[0] * b[2]) * c[1] +
                              (a[0] * b[1] - a[1] * b[0]) * c[2];
        if (triple >= 0)
        {
            continue;
        }

        auto* const node1 = e->getNode(1);
        auto* const node4 = e->getNode(4);
        e->setNode(1, e->getNode(2));
        e->setNode(2, node1);
        e->setNode(4, e->getNode(5));
        e->setNode(5, node4);
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

    fixPrismNodeOrder(*ogs_grid);

    checkCellTypes(*ogs_grid, cell_type);
    checkVolume(*ogs_grid);

    return ogs_grid;
}
