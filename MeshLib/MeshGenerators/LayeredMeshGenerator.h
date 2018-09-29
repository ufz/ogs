/**
 * \file   LayeredMeshGenerator.h
 * \author Karsten Rink
 * \date   2014-09-18
 * \brief  Definition of the SubsurfaceMapper class
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <memory>
#include <limits>
#include <string>
#include <vector>

namespace GeoLib {
    class Raster;
}

namespace MeshLib {
    class Mesh;
    class Node;
    class Element;
}

/**
 * \brief Base class for creation of 3D subsurface meshes based on raster data
 */
class LayeredMeshGenerator
{
public:
    /**
    * Returns a subsurface representation of a region represented by a 2D  mesh by reading raster files and calling the appropriate construction method.
    * @param mesh                    The 2D surface mesh that is used as a basis for the subsurface mesh
    * @param rasters                 Containing all the rasters for the subsurface layers from bottom to top (starting with the bottom of the oldest layer and ending with the DEM)
    * @param minimum_thickness       Minimum thickness of each of the newly created layers (i.e. nodes with a vertical distance smaller than this will be collapsed)
    * @param noDataReplacementValue  Default z-coordinate if there are mesh nodes not located on the DEM raster (i.e. raster_paths[0])
    * @result true if the subsurface representation has been created, false if there was an error
    */
    virtual bool createLayers(MeshLib::Mesh const& mesh,
                              std::vector<GeoLib::Raster const*> const& rasters,
                              double minimum_thickness,
                              double noDataReplacementValue = 0.0) final;

    /**
    * Constructs a subsurface representation based on a 2D mesh and a number of rasters representing subsurface layer boundaries.
    * @param mesh                    The 2D surface mesh that is used as a basis for the subsurface mesh
    * @param rasters                 Containing all the raster-data for the subsurface layers from bottom to top (starting with the bottom of the oldest layer and ending with the DEM)
    * @param minimum_thickness       Minimum thickness of each of the newly created layers (i.e. nodes with a vertical distance smaller than this will be collapsed)
    * @param noDataReplacementValue  Default z-coordinate if there are mesh nodes not located on the DEM raster (i.e. raster_paths[0])
    * @result true if the subsurface representation has been created, false if there was an error
    */
    virtual bool createRasterLayers(MeshLib::Mesh const& mesh,
                                    std::vector<GeoLib::Raster const*> const& rasters,
                                    double minimum_thickness,
                                    double noDataReplacementValue) = 0;

    /// Returns a mesh of the subsurface representation
    std::unique_ptr<MeshLib::Mesh> getMesh(std::string const& mesh_name) const;

protected:
    LayeredMeshGenerator();
    virtual ~LayeredMeshGenerator() = default;

    /// Adds another layer to the subsurface mesh
    virtual void addLayerToMesh(MeshLib::Mesh const& mesh_layer, unsigned layer_id, GeoLib::Raster const& raster) = 0;

    /**
    * Calculates the node Position of a subsurface node based on the given raster but also constrained by the DEM layer
    * as an upper bound and the the layer located below as a lower bound (i.e. older stratigraphic layers are favored and
    * nodes cannot be located above surface).
    * @param dem_node          The node at this xy-location on the DEM
    * @param last_layer_node   The node at this xy-location on the layer below
    * @param raster            The raster file for the current layer
    * @param new_node_id       Node ID to be used if there is a meaningful node position to be found
    * @result A new node at the given xy position.
    */
    MeshLib::Node* getNewLayerNode(MeshLib::Node const& dem_node,
                                   MeshLib::Node const& last_layer_node,
                                   GeoLib::Raster const& raster,
                                   std::size_t new_node_id) const;

    /// Calculates a data-dependent epsilon value
    double calcEpsilon(GeoLib::Raster const& low, GeoLib::Raster const& high);

    /// Cleans up the already created objects in case of an error
    void cleanUpOnError();

    double _elevation_epsilon;
    double _minimum_thickness;
    std::vector<int> _materials;
    std::vector<MeshLib::Node*> _nodes;
    std::vector<MeshLib::Element*> _elements;
};
