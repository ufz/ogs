/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "RasterToMesh.h"

#include <vtkCell.h>
#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>

#include "BaseLib/Logging.h"
#include "MeshGenerator.h"
#include "MeshLib/Elements/Elements.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshSearch/ElementSearch.h"
#include "MeshLib/Node.h"
#include "MeshToolsLib/MeshEditing/RemoveMeshComponents.h"
#include "MeshToolsLib/MeshGenerators/MeshGenerator.h"

namespace MeshToolsLib
{
std::unique_ptr<MeshLib::Mesh> RasterToMesh::convert(
    GeoLib::Raster const& raster,
    MeshLib::MeshElemType elem_type,
    MeshLib::UseIntensityAs intensity_type,
    std::string const& array_name)
{
    return convert(raster.data(),
                   raster.getHeader(),
                   elem_type,
                   intensity_type,
                   array_name);
}

std::unique_ptr<MeshLib::Mesh> RasterToMesh::convert(
    double const* const img,
    GeoLib::RasterHeader const& header,
    MeshLib::MeshElemType elem_type,
    MeshLib::UseIntensityAs intensity_type,
    std::string const& array_name)
{
    if ((elem_type != MeshLib::MeshElemType::TRIANGLE) &&
        (elem_type != MeshLib::MeshElemType::QUAD) &&
        (elem_type != MeshLib::MeshElemType::HEXAHEDRON) &&
        (elem_type != MeshLib::MeshElemType::PRISM))
    {
        ERR("Invalid Mesh Element Type.");
        return nullptr;
    }

    if (((elem_type == MeshLib::MeshElemType::TRIANGLE) ||
         (elem_type == MeshLib::MeshElemType::QUAD)) &&
        header.n_depth != 1)
    {
        ERR("Triangle or Quad elements cannot be used to construct meshes from "
            "3D rasters.");
        return nullptr;
    }

    if (intensity_type == MeshLib::UseIntensityAs::ELEVATION &&
        ((elem_type == MeshLib::MeshElemType::PRISM) ||
         (elem_type == MeshLib::MeshElemType::HEXAHEDRON)))
    {
        ERR("Elevation mapping can only be performed for 2D meshes.");
        return nullptr;
    }

    std::unique_ptr<MeshLib::Mesh> mesh;
    if (elem_type == MeshLib::MeshElemType::TRIANGLE)
    {
        mesh.reset(MeshToolsLib::MeshGenerator::generateRegularTriMesh(
            header.n_cols,
            header.n_rows,
            header.cell_size,
            header.origin,
            "RasterDataMesh"));
    }
    else if (elem_type == MeshLib::MeshElemType::QUAD)
    {
        mesh.reset(MeshToolsLib::MeshGenerator::generateRegularQuadMesh(
            header.n_cols,
            header.n_rows,
            header.cell_size,
            header.origin,
            "RasterDataMesh"));
    }
    else if (elem_type == MeshLib::MeshElemType::PRISM)
    {
        mesh.reset(MeshToolsLib::MeshGenerator::generateRegularPrismMesh(
            header.n_cols,
            header.n_rows,
            header.n_depth,
            header.cell_size,
            header.origin,
            "RasterDataMesh"));
    }
    else if (elem_type == MeshLib::MeshElemType::HEXAHEDRON)
    {
        mesh.reset(MeshToolsLib::MeshGenerator::generateRegularHexMesh(
            header.n_cols,
            header.n_rows,
            header.n_depth,
            header.cell_size,
            header.origin,
            "RasterDataMesh"));
    }

    std::unique_ptr<MeshLib::Mesh> new_mesh;
    std::vector<std::size_t> elements_to_remove;
    if (intensity_type == MeshLib::UseIntensityAs::ELEVATION)
    {
        std::vector<MeshLib::Node*> const& nodes(mesh->getNodes());
        std::vector<MeshLib::Element*> const& elems(mesh->getElements());
        std::size_t const n_nodes(elems[0]->getNumberOfNodes());
        bool const double_idx =
            (elem_type == MeshLib::MeshElemType::TRIANGLE) ||
            (elem_type == MeshLib::MeshElemType::PRISM);
        std::size_t const m = (double_idx) ? 2 : 1;
        for (std::size_t k = 0; k < header.n_depth; k++)
        {
            std::size_t const layer_idx = (k * header.n_rows * header.n_cols);
            for (std::size_t i = 0; i < header.n_cols; i++)
            {
                std::size_t const idx(i * header.n_rows + layer_idx);
                for (std::size_t j = 0; j < header.n_rows; j++)
                {
                    double const val(img[idx + j]);
                    if (val == header.no_data)
                    {
                        elements_to_remove.push_back(m * (idx + j));
                        if (double_idx)
                        {
                            elements_to_remove.push_back(m * (idx + j) + 1);
                        }
                        continue;
                    }
                    for (std::size_t n = 0; n < n_nodes; ++n)
                    {
                        (*(nodes[getNodeIndex(*elems[m * (idx + j)], n)]))[2] =
                            val;
                        if (double_idx)
                        {
                            (*(nodes[getNodeIndex(*elems[m * (idx + j) + 1],
                                                  n)]))[2] = val;
                        }
                    }
                }
            }
        }
    }
    else
    {
        MeshLib::Properties& properties = mesh->getProperties();
        MeshLib::ElementSearch ex(*mesh);
        if (array_name == "MaterialIDs")
        {
            auto* const prop_vec = properties.createNewPropertyVector<int>(
                array_name, MeshLib::MeshItemType::Cell, 1);
            fillPropertyVector<int>(*prop_vec, img, header, elem_type);
            ex.searchByPropertyValue<int>(array_name,
                                          static_cast<int>(header.no_data));
        }
        else
        {
            auto* const prop_vec = properties.createNewPropertyVector<double>(
                array_name, MeshLib::MeshItemType::Cell, 1);
            fillPropertyVector<double>(*prop_vec, img, header, elem_type);
            ex.searchByPropertyValue<double>(array_name, header.no_data);
        }
        elements_to_remove = ex.getSearchedElementIDs();
    }
    if (!elements_to_remove.empty())
    {
        new_mesh.reset(MeshToolsLib::removeElements(
            *mesh, elements_to_remove, mesh->getName()));
    }
    else
    {
        new_mesh = std::move(mesh);
    }

    if (intensity_type == MeshLib::UseIntensityAs::NONE)
    {
        new_mesh->getProperties().removePropertyVector(array_name);
    }

    return new_mesh;
}

}  // namespace MeshToolsLib
