/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 */

#pragma once

#include <logog/include/logog.hpp>

#include "GeoLib/Raster.h"
#include "MeshLib/Location.h"
#include "MeshLib/MeshEnums.h"
#include "MeshLib/Properties.h"

class vtkImageData; // For conversion from Image to QuadMesh

namespace MeshLib {

class Mesh;
class Node;
class Element;
template <typename T>
class PropertyVector;

/**
 * \brief Converts raster data into an OGS mesh.
 */
class RasterToMesh
{
public:
    /**
     * Converts greyscale raster into a mesh.
     * \param raster         input raster.
     * \param elem_type      defines if elements of the new mesh should be
     *                       triangles or quads (or hexes for 3D).
     * \param intensity_type defines how image intensities are interpreted.
     * \param array_name     mesh property name, defaults to "Colour" if not
     *                       given.
     */
    static MeshLib::Mesh* convert(GeoLib::Raster const& raster,
                                  MeshElemType elem_type,
                                  UseIntensityAs intensity_type,
                                  std::string const& array_name = "Colour");

    /**
     * Converts a vtkImageData into a mesh.
     * \param img            input image.
     * \param origin         coordinates of image's origin, lower left corner.
     * \param scalingFactor  edge length of each pixel
     * \param elem_type      defines if elements of the new mesh should be
     *                       triangles or quads (or hexes for 3D).
     * \param intensity_type defines how image intensities are interpreted.
     * \param array_name     mesh property name, defaults to "Colour" if not
     *                       given.
     */
    static MeshLib::Mesh* convert(vtkImageData* img,
                                  const double origin[3],
                                  const double scalingFactor,
                                  MeshElemType elem_type,
                                  UseIntensityAs intensity_type,
                                  std::string const& array_name = "Colour");

    /**
     * Converts double array with raster values into a mesh.
     * \param img            input image.
     * \param header         raster header information.
     * \param elem_type      defines if elements of the new mesh should be
     *                       triangles or quads (or hexes for 3D).
     * \param intensity_type defines how image intensities are interpreted.
     * \param array_name     mesh property name, defaults to "Colour" if not
     *                       given.
     */
    static MeshLib::Mesh* convert(const double* img,
                                  GeoLib::RasterHeader const& header,
                                  MeshElemType elem_type,
                                  UseIntensityAs intensity_type,
                                  std::string const& array_name = "Colour");

private:
    static MeshLib::Mesh* constructMesh(
        std::vector<double> const& pix_val,
        std::string const& array_name,
        std::vector<bool> const& pix_vis,
        GeoLib::RasterHeader const& header,
        MeshLib::MeshElemType elem_type,
        MeshLib::UseIntensityAs intensity_type);

    static std::vector<MeshLib::Node*> createNodeVector(
        std::vector<double> const& elevation,
        std::vector<int> &node_idx_map,
        GeoLib::RasterHeader const& header,
        bool use_elevation);

    static std::vector<MeshLib::Element*> createElementVector(
        std::vector<bool> const& pix_vis,
        std::vector<MeshLib::Node*> const& nodes,
        std::vector<int> const& node_idx_map,
        std::size_t const imgHeight,
        std::size_t const imgWidth,
        MeshElemType elem_type);

    template<typename T>
    static void fillPropertyVector(
        MeshLib::PropertyVector<T> &prop_vec,
        std::vector<double> const& pix_val,
        std::vector<bool> const& pix_vis,
        const std::size_t &imgHeight,
        const std::size_t &imgWidth,
        MeshElemType elem_type)
    {
        for (std::size_t i = 0; i < imgHeight; i++)
            for (std::size_t j = 0; j < imgWidth; j++)
            {
                if (!pix_vis[i*imgWidth+j])
                    continue;
                T val (static_cast<T>(pix_val[i*(imgWidth+1)+j]));
                if (elem_type == MeshElemType::TRIANGLE)
                {
                    prop_vec.push_back(val);
                    prop_vec.push_back(val);
                }
                else if (elem_type == MeshElemType::QUAD)
                    prop_vec.push_back(val);
            }
    }

    static double getExistingValue(const double* img, std::size_t length);
};

} // end namespace MeshLib
