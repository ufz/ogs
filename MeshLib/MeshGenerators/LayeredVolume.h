/**
 * \file   LayeredVolume.h
 * \author Karsten Rink
 * \date   2014-04-11
 * \brief  Definition of the LayeredVolume class
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <string>
#include <vector>

#include "MeshLib/Node.h"
#include "LayeredMeshGenerator.h"

namespace GeoLib {
    class GEOObjects;
    class Surface;
}

/**
 * \brief Creates a volume geometry from 2D mesh layers based on raster data.
 */
class LayeredVolume : public LayeredMeshGenerator
{
public:
    LayeredVolume() {}
    ~LayeredVolume() {}

    /**
     * Constructs a subsurface representation of a mesh using only 2D elements (i.e. layer boundaries are represented by surfaces)
     * @param mesh                    The 2D surface mesh that is used as a basis for the subsurface mesh
     * @param rasters                 Containing all the raster-data for the subsurface layers from bottom to top (starting with the bottom of the oldest layer and ending with the DEM)
     * @param minimum_thickness       Minimum thickness of each of the newly created layers (i.e. nodes with a vertical distance smaller than this will be collapsed)
     * @param noDataReplacementValue  Default z-coordinate if there are mesh nodes not located on the DEM raster (i.e. raster_paths[0])
     * @result true if the subsurface representation has been created, false if there was an error
     */
    bool createRasterLayers(const MeshLib::Mesh &mesh,
                            const std::vector<GeoLib::Raster const*> &rasters,
                            double minimum_thickness,
                            double noDataReplacementValue = 0.0);

    /// Returns the region attribute vector necessary for assigning region attributes via TetGen
    std::vector<MeshLib::Node> getAttributePoints() { return _attribute_points; }

private:
    /// Adds another layer to the subsurface mesh
    void addLayerToMesh(const MeshLib::Mesh &mesh_layer, unsigned layer_id, GeoLib::Raster const& raster);

    /// Creates boundary surfaces between the mapped layers to make the volumes watertight
    void addLayerBoundaries(const MeshLib::Mesh &layer, std::size_t nLayers);

    /// Removes duplicate 2D elements (possible due to outcroppings)
    void removeCongruentElements(std::size_t nLayers, std::size_t nElementsPerLayer);

    std::vector<MeshLib::Node> _attribute_points;
};
