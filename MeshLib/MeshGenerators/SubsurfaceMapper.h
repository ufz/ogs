/**
 * \file   SubsurfaceMapper.h
 * \author Karsten Rink
 * \date   2014-09-18
 * \brief  Definition of the SubsurfaceMapper class
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef SUBSURFACEMAPPER_H
#define SUBSURFACEMAPPER_H

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
class SubsurfaceMapper
{
public:
    /**
     * Returns a subsurface representation of a region represented by a 2D by reading raster files and calling the appropriate construction method.
     * @param mesh                    The 2D surface mesh that is used as a basis for the subsurface mesh
     * @param raster_paths            Containing all the raster-file-names for the subsurface layers from bottom to top (starting with the bottom of the oldest layer and ending with the DEM)
     * @param noDataReplacementValue  Default z-coordinate if there are mesh nodes not located on the DEM raster (i.e. raster_paths[0]) 
     * @result true if the subsurface representation has been created, false if there was an error
     */
    virtual bool createLayers(MeshLib::Mesh const& mesh, std::vector<std::string> const& raster_paths, double noDataReplacementValue = 0.0) final;

    /**
     * Constructs a subsurface representation based on a 2D mesh and a number of rasters representing subsurface layer boundaries.
     * @param mesh                    The 2D surface mesh that is used as a basis for the subsurface mesh
     * @param rasters                 Containing all the raster-data for the subsurface layers from bottom to top (starting with the bottom of the oldest layer and ending with the DEM)
     * @param noDataReplacementValue  Default z-coordinate if there are mesh nodes not located on the DEM raster (i.e. raster_paths[0]) 
     * @result true if the subsurface representation has been created, false if there was an error
     */
    virtual bool createRasterLayers(MeshLib::Mesh const& mesh, std::vector<GeoLib::Raster const*> const& rasters, double noDataReplacementValue) = 0;

    /// Returns a mesh of the subsurface representation
    MeshLib::Mesh* getMesh(std::string const& mesh_name) const;

protected:
    SubsurfaceMapper();
    ~SubsurfaceMapper() {}

    /// Adds another layer to the subsurface mesh
    virtual void addLayerToMesh(MeshLib::Mesh const& mesh_layer, unsigned layer_id, GeoLib::Raster const& raster) = 0;

    /// Calculates a data-dependent epsilon value
    double calcEpsilon(GeoLib::Raster const& high, GeoLib::Raster const& low);

    /// Checks if all raster files actually exist
    bool allRastersExist(std::vector<std::string> const& raster_paths) const;

    /// Cleans up the already created objects in case of an error
    void cleanUpOnError();

    double _elevation_epsilon;
    std::vector<MeshLib::Node*> _nodes;
    std::vector<MeshLib::Element*> _elements;
};

#endif //SUBSURFACEMAPPER_H
