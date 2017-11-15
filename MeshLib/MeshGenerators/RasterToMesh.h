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
    static MeshLib::Mesh* convert(const double*const img,
        GeoLib::RasterHeader const& header,
        MeshElemType elem_type,
        UseIntensityAs intensity_type,
        std::string const& array_name = "Colour");

private:
    template<typename T>
    static void fillPropertyVector(
        MeshLib::PropertyVector<T> &prop_vec,
        double const*const img,
        GeoLib::RasterHeader const& header,
        MeshElemType elem_type)
    {
        for (std::size_t i = 0; i < header.n_cols; i++)
        {
            std::size_t idx(i * header.n_rows);
            for (std::size_t j = 0; j < header.n_rows; j++)
            {
                auto val(static_cast<T>(img[idx + j]));
                prop_vec.push_back(val);
                if (elem_type == MeshElemType::TRIANGLE || elem_type == MeshElemType::PRISM)
                    prop_vec.push_back(val); // because each pixel is represented by two cells
            }
        }
    }
};

} // end namespace MeshLib
