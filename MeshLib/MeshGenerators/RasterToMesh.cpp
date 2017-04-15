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
        {static_cast<std::size_t>(dims[0]), static_cast<std::size_t>(dims[1]), orig, scalingFactor, -9999};
    const std::size_t incHeight = header.n_rows+1;
    const std::size_t incWidth  = header.n_cols+1;
    std::vector<double> pix_val (incHeight * incWidth, std::numeric_limits<double>::max());
    std::vector<bool> pix_vis (header.n_rows * header.n_cols, false);

    for (std::size_t i = 0; i < header.n_rows; i++)
        for (std::size_t j = 0; j < header.n_cols; j++)
        {
            std::size_t const img_idx = i*header.n_cols + j;
            std::size_t const fld_idx = i*incWidth + j;

            // colour of current pixel
            double* colour = pixelData->GetTuple(img_idx);
            // is current pixel visible?
            bool const visible = (nTuple == 2 || nTuple == 4) ? (colour[nTuple-1] != 0) : true;
            if (!visible)
                continue;

            double const value = (nTuple < 3) ?
                colour[0] : // grey (+ alpha)
                (0.3 * colour[0] + 0.6 * colour[1] + 0.1 * colour[2]); // rgb(a)
            pix_vis[img_idx] = true;
            pix_val[fld_idx] = value;
            pix_val[fld_idx+1] = value;
            pix_val[fld_idx+incWidth] = value;
            pix_val[fld_idx+incWidth+1] = value;
        }

    return constructMesh(pix_val, array_name, pix_vis, header, elem_type, intensity_type);
}

MeshLib::Mesh* RasterToMesh::convert(
    double const* img,
    GeoLib::RasterHeader const& header,
    MeshElemType elem_type,
    UseIntensityAs intensity_type,
    std::string const& array_name)
{
    if ((elem_type != MeshElemType::TRIANGLE) && (elem_type != MeshElemType::QUAD))
    {
        ERR("Invalid Mesh Element Type.");
        return nullptr;
    }

    std::size_t const incHeight (header.n_rows+1);
    std::size_t const incWidth (header.n_cols+1);
    std::vector<double> pix_val (incHeight * incWidth, std::numeric_limits<double>::max());
    std::vector<bool> pix_vis (incHeight * incWidth, false);

    for (std::size_t i = 0; i < header.n_rows; i++)
        for (std::size_t j = 0; j < header.n_cols; j++)
        {
            std::size_t const img_idx = i*header.n_cols + j;
            std::size_t const fld_idx = i*incWidth + j;
            if (img[img_idx] == -9999)
                continue;

            pix_vis[img_idx] = true;
            pix_val[fld_idx] = img[img_idx];
            pix_val[fld_idx+1] = img[img_idx];
            pix_val[fld_idx+incWidth] = img[img_idx];
            pix_val[fld_idx+incWidth+1] = img[img_idx];
        }

    return constructMesh(pix_val, array_name, pix_vis, header, elem_type, intensity_type);
}

MeshLib::Mesh* RasterToMesh::constructMesh(
    std::vector<double> const& pix_val,
    std::string const& array_name,
    std::vector<bool> const& pix_vis,
    GeoLib::RasterHeader const& header,
    MeshLib::MeshElemType elem_type,
    MeshLib::UseIntensityAs intensity_type)
{
    std::vector<int> node_idx_map ((header.n_rows+1) * (header.n_cols+1), -1);
    bool const use_elevation (intensity_type == MeshLib::UseIntensityAs::ELEVATION);
    std::vector<MeshLib::Node*> nodes (createNodeVector(pix_val, node_idx_map, header, use_elevation));
    if (nodes.empty())
        return nullptr;

    std::vector<MeshLib::Element*> elements(createElementVector(
        pix_vis, nodes, node_idx_map, header.n_rows, header.n_cols, elem_type));
    if (elements.empty())
        return nullptr;

    MeshLib::Properties properties;
    if (intensity_type == MeshLib::UseIntensityAs::MATERIALS)
    {
        auto* const prop_vec = properties.createNewPropertyVector<int>(
            "MaterialIDs", MeshLib::MeshItemType::Cell, 1);
        fillPropertyVector<int>(*prop_vec, pix_val, pix_vis, header.n_rows, header.n_cols, elem_type);
    }
    else if (intensity_type == MeshLib::UseIntensityAs::DATAVECTOR)
    {
        auto* const prop_vec = properties.createNewPropertyVector<double>(
            array_name, MeshLib::MeshItemType::Cell, 1);
        fillPropertyVector<double>(*prop_vec, pix_val, pix_vis, header.n_rows, header.n_cols, elem_type);
    }

    return new MeshLib::Mesh("RasterDataMesh", nodes, elements, properties);
}

std::vector<MeshLib::Node*> RasterToMesh::createNodeVector(
    std::vector<double> const&  elevation,
    std::vector<int> & node_idx_map,
    GeoLib::RasterHeader const& header,
    bool use_elevation)
{
    std::size_t node_idx_count(0);
    double const x_offset(header.origin[0] - header.cell_size/2.0);
    double const y_offset(header.origin[1] - header.cell_size/2.0);
    std::vector<MeshLib::Node*> nodes;
    for (std::size_t i = 0; i < (header.n_rows+1); i++)
        for (std::size_t j = 0; j < (header.n_cols+1); j++)
        {
            std::size_t const index = i * (header.n_cols+1) + j;
            if (elevation[index] == std::numeric_limits<double>::max())
                continue;

            double const zValue = (use_elevation) ? elevation[index] : 0;
            auto* node(new MeshLib::Node(x_offset + (header.cell_size * j),
                                         y_offset + (header.cell_size * i),
                                         zValue));
            nodes.push_back(node);
            node_idx_map[index] = node_idx_count;
            node_idx_count++;
        }
    return nodes;
}

std::vector<MeshLib::Element*> RasterToMesh::createElementVector(
    std::vector<bool> const& pix_vis,
    std::vector<MeshLib::Node*> const& nodes,
    std::vector<int> const&node_idx_map,
    std::size_t const imgHeight,
    std::size_t const imgWidth,
    MeshElemType elem_type)
{
    std::vector<MeshLib::Element*> elements;
    std::size_t const incWidth (imgWidth+1);
    for (std::size_t i = 0; i < imgHeight; i++)
        for (std::size_t j = 0; j < imgWidth; j++)
        {
            if (!pix_vis[i*imgWidth+j])
                continue;

            int const idx = i * incWidth + j;
            if (elem_type == MeshElemType::TRIANGLE)
            {
                auto** tri1_nodes = new MeshLib::Node*[3];
                tri1_nodes[0] = nodes[node_idx_map[idx]];
                tri1_nodes[1] = nodes[node_idx_map[idx+1]];
                tri1_nodes[2] = nodes[node_idx_map[idx+incWidth]];

                auto** tri2_nodes = new MeshLib::Node*[3];
                tri2_nodes[0] = nodes[node_idx_map[idx+1]];
                tri2_nodes[1] = nodes[node_idx_map[idx+incWidth+1]];
                tri2_nodes[2] = nodes[node_idx_map[idx+incWidth]];

                elements.push_back(new MeshLib::Tri(tri1_nodes)); // upper left triangle
                elements.push_back(new MeshLib::Tri(tri2_nodes)); // lower right triangle
            }
            else if (elem_type == MeshElemType::QUAD)
            {
                auto** quad_nodes = new MeshLib::Node*[4];
                quad_nodes[0] = nodes[node_idx_map[idx]];
                quad_nodes[1] = nodes[node_idx_map[idx + 1]];
                quad_nodes[2] = nodes[node_idx_map[idx + incWidth + 1]];
                quad_nodes[3] = nodes[node_idx_map[idx + incWidth]];
                elements.push_back(new MeshLib::Quad(quad_nodes));
            }
        }
    return elements;
}

double RasterToMesh::getExistingValue(const double* img, std::size_t length)
{
    for (std::size_t i=0; i<length; i++)
    {
        if (img[i] != -9999)
            return img[i];
    }
    return -9999;
}

} // end namespace MeshLib
