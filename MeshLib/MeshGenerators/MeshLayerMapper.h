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
    * Based on a triangle-or quad mesh this method creates a 3D mesh with with a given number of prism- or hex-layers
    * \param mesh The triangle/quad mesh that is the basis for the new prism/hex mesh
    * \param nLayers The number of layers of prism/hex elements that will be extruded from the triangle/quad elements of the original mesh
    * \param thickness The thickness of each of these newly added layers
    * \param mesh_name The name of the newly created mesh
    * \return A mesh with the requested number of layers of prism/hex elements
    */
    MeshLib::Mesh* createStaticLayers(MeshLib::Mesh const& mesh, std::vector<float> const& layer_thickness_vector, std::string const& mesh_name = "SubsurfaceMesh") const;

    MeshLib::Mesh* createRasterLayers(MeshLib::Mesh const& mesh, std::vector<std::string> const& raster_paths, std::string const& mesh_name = "SubsurfaceMesh") const;

    MeshLib::Mesh* createRasterLayers(MeshLib::Mesh const& mesh, std::vector<GeoLib::Raster const*> const& rasters, std::string const& mesh_name = "SubsurfaceMesh") const;

    /**
    * Maps the z-values of nodes in the designated layer of the given mesh according to the given raster.
    * Note: This only results in a valid mesh if the layers don't intersect each other.
    */
    bool layerMapping(MeshLib::Mesh &mesh, const std::string &rasterfile,
                      const unsigned nLayers, const unsigned layer_id, double noDataReplacementValue) const;

    /**
    * Maps the z-values of nodes in the designated layer of the given mesh according to the given raster.
    * Note: This only results in a valid mesh if the layers don't intersect each other.
    */
    bool layerMapping(MeshLib::Mesh &mesh, const GeoLib::Raster &raster,
                      const unsigned nLayers, const unsigned layer_id, double noDataReplacementValue) const;

    /**
    * Blends a mesh with the surface given by dem_raster. Nodes and elements above the surface are either removed or adapted to fit the surface.
    * Note: It is unlikely but possible that the new nodes vector contains (very few) nodes that are not part of any element. This problem is
    * remedied at the end of method upon creating the actual mesh from the new node- and element-vector as the mesh-constructor checks for such
    * nodes and removes them. This note is just to call this issue to attention in case this methods is changed.
    */
    MeshLib::Mesh* blendLayersWithSurface(MeshLib::Mesh &mesh, const unsigned nLayers, const std::string &dem_raster) const;

private:
    /// Adds another layer to the subsurface mesh
    void addLayerToMesh(const MeshLib::Mesh &mesh_layer, unsigned layer_id, GeoLib::Raster const& raster, std::vector<MeshLib::Node*> &new_nodes, std::vector<MeshLib::Element*> &new_elems) const;

    /// Checks if all raster files actually exist
	bool allRastersExist(const std::vector<std::string> &raster_paths) const;

    static const unsigned _pyramid_base[3][4];
};

#endif //MESHLAYERMAPPER_H
