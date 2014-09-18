/**
 * \file   MeshLayerMapper.h
 * \author Karsten Rink
 * \date   2010-11-01
 * \brief  Definition of the MeshLayerMapper class.
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MESHLAYERMAPPER_H
#define MESHLAYERMAPPER_H

#include <string>
#include "GeoLib/Raster.h"

class QImage;

namespace MeshLib {
    class Mesh;
    class Node;
    class Element;
}

/**
* \brief Manipulating and adding layers to an existing mesh
*/
class MeshLayerMapper
{
public:
    MeshLayerMapper() {};
    ~MeshLayerMapper() {};

    /**
    * Based on a 2D triangle-or quad mesh this method creates a 3D mesh with with a given number of prism- or hex-layers
    * \param mesh The triangle/quad mesh that is the basis for the new prism/hex mesh
    * \param nLayers The number of layers of prism/hex elements that will be extruded from the triangle/quad elements of the original mesh
    * \param thickness The thickness of each of these newly added layers
    * \param mesh_name The name of the newly created mesh
    * \return A mesh with the requested number of layers of prism/hex elements
    */
    MeshLib::Mesh* createStaticLayers(MeshLib::Mesh const& mesh, std::vector<float> const& layer_thickness_vector, std::string const& mesh_name = "SubsurfaceMesh") const;

   /**
    * Based on a 2D triangle mesh this method creates a 3D mesh with with a given number of prism-layers.
    * Note: This method would theoretically also work with quad meshes. However, this is discouraged as quad elements will most likely not
    * be coplanar after the mapping process and result in invaled elements.
    * \param mesh The 2D triangle mesh that is the basis for the new 3D prism mesh
    * \param nLayers The number of layers of prism/hex elements that will be extruded from the triangle/quad elements of the original mesh
    * \param thickness The thickness of each of these newly added layers
    * \param mesh_name The name of the newly created mesh
    * \return A mesh with the requested number of layers of prism/hex elements
    */
    MeshLib::Mesh* createRasterLayers(MeshLib::Mesh const& mesh, std::vector<std::string> const& raster_paths, std::string const& mesh_name = "SubsurfaceMesh") const;

    MeshLib::Mesh* createRasterLayers(MeshLib::Mesh const& mesh, std::vector<GeoLib::Raster const*> const& rasters, std::string const& mesh_name = "SubsurfaceMesh") const;

    /**
    * Maps the elevation of nodes of a given 2D mesh according to the raster specified by the file path.
    * At locations wher no information is given, node elevation is set to noDataReplacementValue.
    */
    bool layerMapping(MeshLib::Mesh &mesh, const std::string &rasterfile, double noDataReplacementValue) const;

    /**
    * Maps the elevation of nodes of a given 2D mesh according to the raster. At locations wher no 
    * information is given, node elevation is set to noDataReplacementValue.
    */
    bool layerMapping(MeshLib::Mesh &mesh, const GeoLib::Raster &raster, double noDataReplacementValue) const;

private:
    /// Adds another layer to a subsurface mesh
    void addLayerToMesh(const MeshLib::Mesh &mesh_layer, unsigned layer_id, GeoLib::Raster const& raster, std::vector<MeshLib::Node*> &new_nodes, std::vector<MeshLib::Element*> &new_elems) const;

    /// Checks if all raster files actually exist
	bool allRastersExist(const std::vector<std::string> &raster_paths) const;

    static const unsigned _pyramid_base[3][4];
};

#endif //MESHLAYERMAPPER_H
