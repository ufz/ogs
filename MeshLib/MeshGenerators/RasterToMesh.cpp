/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "RasterToMesh.h"

#include "MeshLib/Elements/Elements.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"
#include "MeshLib/MeshEditing/RemoveMeshComponents.h"
#include "MeshLib/MeshGenerators/MeshGenerator.h"
#include "MeshLib/MeshSearch/ElementSearch.h"

#include <vtkImageData.h>
#include <vtkDataArray.h>
#include <vtkCell.h>
#include <vtkCellData.h>
#include <vtkSmartPointer.h>
#include <vtkPointData.h>

namespace MeshLib
{

MeshLib::Mesh* RasterToMesh::convert(
    GeoLib::Raster const& raster,
    MeshElemType elem_type,
    UseIntensityAs intensity_type,
    std::string const& array_name)
{
    //return convert(raster.begin(), raster.getHeader(), elem_type, intensity_type, array_name);
    return convert(raster.begin(), raster.getHeader(), elem_type, intensity_type, array_name);
}

MeshLib::Mesh* RasterToMesh::convert(
    vtkImageData* img,
    const double origin[3],
    const double scalingFactor,
    MeshElemType elem_type,
    UseIntensityAs intensity_type,
    std::string const& array_name)
{
    if ((elem_type != MeshElemType::TRIANGLE) && (elem_type != MeshElemType::QUAD))
    {
        ERR("VtkMeshConverter::convertImgToMesh(): Invalid Mesh Element Type.");
        return nullptr;
    }

    vtkSmartPointer<vtkDataArray> pixelData = vtkSmartPointer<vtkDataArray>(img->GetPointData()->GetScalars());
    int* dims = img->GetDimensions();
    int nTuple = pixelData->GetNumberOfComponents();
    if (nTuple < 1 || nTuple > 4)
    {
        ERR("VtkMeshConverter::convertImgToMesh(): Unsupported pixel composition!");
        return nullptr;
    }

    MathLib::Point3d const orig (std::array<double,3>{{origin[0], origin[1], origin[2]}});
    GeoLib::RasterHeader const header =
        {static_cast<std::size_t>(dims[0]), static_cast<std::size_t>(dims[1]), 1, orig, scalingFactor, -9999};

    double *pix = new double [header.n_cols * header.n_rows];
    for (std::size_t i = 0; i < header.n_rows; i++)
    {
        std::size_t const idx = i*header.n_cols;
        for (std::size_t j = 0; j < header.n_cols; j++)
        {
            double* colour = pixelData->GetTuple(idx+j);
            bool const visible = (nTuple == 2 || nTuple == 4) ? (colour[nTuple - 1] != 0) : true;
            if (!visible)
                pix[idx+j] = header.no_data;
            else 
                pix[idx + j] = (nTuple < 3) ?
                    colour[0] : // grey (+ alpha)
                    (0.3 * colour[0] + 0.6 * colour[1] + 0.1 * colour[2]); // rgb(a)
        }
    }

    return convert(pix, header, elem_type, intensity_type, array_name);
}

MeshLib::Mesh* RasterToMesh::convert(
    double const*const img,
    GeoLib::RasterHeader const& header,
    MeshElemType elem_type,
    UseIntensityAs intensity_type,
    std::string const& array_name)
{
    if ((elem_type != MeshElemType::TRIANGLE) && 
        (elem_type != MeshElemType::QUAD) &&
        (elem_type != MeshElemType::HEXAHEDRON) &&
        (elem_type != MeshElemType::PRISM))
    {
        ERR("Invalid Mesh Element Type.");
        return nullptr;
    }

    if (((elem_type == MeshElemType::TRIANGLE) || (elem_type == MeshElemType::QUAD)) &&
        header.n_depth != 1)
    {
        ERR("Triangle or Quad elements cannot be used to construct meshes from 3D rasters.");
        return nullptr;
    }

    MathLib::Point3d mesh_origin({header.origin[0]-(header.cell_size/2.0), header.origin[1]-(header.cell_size / 2.0), header.origin[2]});
    std::unique_ptr<MeshLib::Mesh> mesh (nullptr);
    if (elem_type == MeshElemType::TRIANGLE)
        mesh.reset(MeshLib::MeshGenerator::generateRegularTriMesh(header.n_cols, header.n_rows, header.cell_size, mesh_origin, "RasterDataMesh"));
    else if (elem_type == MeshElemType::QUAD)
        mesh.reset(MeshLib::MeshGenerator::generateRegularQuadMesh(header.n_cols, header.n_rows, header.cell_size, mesh_origin, "RasterDataMesh"));
    else if (elem_type == MeshElemType::PRISM)
        mesh.reset(MeshLib::MeshGenerator::generateRegularPrismMesh(header.n_cols, header.n_rows, header.n_depth, header.cell_size, mesh_origin, "RasterDataMesh"));
    else if (elem_type == MeshElemType::HEXAHEDRON)
        mesh.reset(MeshLib::MeshGenerator::generateRegularHexMesh(header.n_cols, header.n_rows, header.n_depth, header.cell_size, mesh_origin, "RasterDataMesh"));

    MeshLib::Mesh* new_mesh(nullptr);
    std::vector<std::size_t> elements_to_remove;
    if (intensity_type == UseIntensityAs::ELEVATION)
    {
        std::vector<MeshLib::Node*> nodes(mesh->getNodes());
        std::vector<MeshLib::Element*> const& elems(mesh->getElements());
        std::size_t const n_nodes(elems[0]->getNumberOfNodes());
        bool const double_idx = (elem_type == MeshElemType::TRIANGLE) || (elem_type == MeshElemType::PRISM);
        std::size_t const m = (double_idx) ? 2 : 1;
        for (std::size_t i = 0; i < header.n_cols; i++)
        {
            std::size_t const idx(i * header.n_rows);
            for (std::size_t j = 0; j < header.n_rows; j++)
            {
                double const val(img[idx + j]);
                if (val == header.no_data)
                {
                    elements_to_remove.push_back(m * (idx + j));
                    if (double_idx)
                        elements_to_remove.push_back(m * (idx + j)+1);
                    continue;
                }
                for (std::size_t n = 0; n < n_nodes; ++n)
                {
                    (*(nodes[elems[m * (idx + j)]->getNodeIndex(n)]))[2] = val;
                    if (double_idx)
                        (*(nodes[elems[m * (idx + j)+1]->getNodeIndex(n)]))[2] = val;
                }
            }
        }
    }
    else
    {
        MeshLib::Properties &properties = mesh->getProperties();
        if (array_name == "MaterialIDs")
        {
            auto* const prop_vec = properties.createNewPropertyVector<int>(
                array_name, MeshLib::MeshItemType::Cell, 1);
            fillPropertyVector<int>(*prop_vec, img, header, elem_type);
        }
        else
        {
            auto* const prop_vec = properties.createNewPropertyVector<double>(
                array_name, MeshLib::MeshItemType::Cell, 1);
        fillPropertyVector<double>(*prop_vec, img, header, elem_type);
        }

        MeshLib::ElementSearch ex(*mesh.get());
        ex.searchByPropertyValue(header.no_data, array_name);
        elements_to_remove = ex.getSearchedElementIDs();
    }
    if (!elements_to_remove.empty())
        new_mesh = MeshLib::removeElements(*mesh.get(), elements_to_remove, mesh->getName());
    else
        new_mesh = mesh.release();

    if (intensity_type == UseIntensityAs::NONE)
        new_mesh->getProperties().removePropertyVector(array_name);

    return new_mesh;
}

} // end namespace MeshLib
