/**
 * \file   LayeredVolume.h
 * \author Karsten Rink
 * \date   2014-04-11
 * \brief  Definition of the LayeredVolume class
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef LAYEREDVOLUME_H
#define LAYEREDVOLUME_H

#include <string>
#include <vector>

#include "Raster.h"
#include "Node.h"

namespace GeoLib {
	class GEOObjects;
	class Surface;
}

namespace MeshLib {
	class Mesh;
	class Element;
}

/**
 * \brief Creates a volume geometry from 2D mesh layers based on raster data.
 */
class LayeredVolume
{
public:
	LayeredVolume();
	~LayeredVolume() {}

	/**
	 * Constructs a subsurface representation of a mesh using only 2D elements (i.e. layer boundaries are represented by surfaces)
	 * @param mesh                    The 2D surface mesh that is used as a basis for the subsurface mesh
	 * @param rasters                 Containing all the raster-data for the subsurface layers from top to bottom (starting with the DEM)
	 * @param noDataReplacementValue  Default z-coordinate if there are mesh nodes not located on the DEM raster (i.e. raster_paths[0]) 
	 * @result true if the subsurface representation has been created, false if there was an error
	 */
	bool createGeoVolumes(const MeshLib::Mesh &mesh, const std::vector<GeoLib::Raster const*> &rasters, double noDataReplacementValue = 0.0);

	/**
	 * Constructs a subsurface representation of a mesh using only 2D elements (i.e. layer boundaries are represented by surfaces)
	 * @param mesh                    The 2D surface mesh that is used as a basis for the subsurface mesh
	 * @param raster_paths            Containing all the raster-file-names for the subsurface layers from top to bottom (starting with the DEM)
	 * @param noDataReplacementValue  Default z-coordinate if there are mesh nodes not located on the DEM raster (i.e. raster_paths[0]) 
	 * @result true if the subsurface representation has been created, false if there was an error
	 */
	bool createGeoVolumes(const MeshLib::Mesh &mesh, const std::vector<std::string> &raster_paths, double noDataReplacementValue = 0.0);

	/// Returns the subsurface representation that can be used as input for the TetGenInterface
	MeshLib::Mesh* getMesh() const { return _mesh; }

	/// Returns the region attribute vector necessary for assigning region attributes via TetGen
	std::vector<MeshLib::Node> getAttributePoints() { return _attribute_points; }

	/// Converts the subsurface representation to geometry objects and adds them to the geo storage object
	bool exportToGeometry(GeoLib::GEOObjects &geo_objects) const;

private:
	/// Adds another layer to the subsurface mesh
	void addLayerToMesh(const MeshLib::Mesh &mesh_layer, unsigned layer_id);

	/// Creates boundary surfaces between the mapped layers to make the volumes watertight
	void addLayerBoundaries(const MeshLib::Mesh &layer, std::size_t nLayers);

	/// Removes duplicate 2D elements (possible due to outcroppings)
	void removeCongruentElements(std::size_t nLayers, std::size_t nElementsPerLayer);

	/// Calculates a data-dependent epsilon value
	double calcEpsilon(const GeoLib::Raster &high, const GeoLib::Raster &low);

	/// Checks if all raster files actually exist
	bool allRastersExist(const std::vector<std::string> &raster_paths) const;

	/// Cleans up the already created object in case of an error
	void cleanUpOnError();

	static const double _invalid_value;
	double _elevation_epsilon;
	std::vector<MeshLib::Node*> _nodes;
	std::vector<MeshLib::Element*> _elements;
	std::vector<MeshLib::Node> _attribute_points;
	MeshLib::Mesh* _mesh;
};

#endif //LAYEREDVOLUME_H
