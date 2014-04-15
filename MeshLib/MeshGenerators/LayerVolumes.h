/**
 * \file   LayerVolumes.h
 * \author Karsten Rink
 * \date   2014-04-11
 * \brief  Definition of the LayerVolumes class
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef LAYERVOLUMES_H
#define LAYERVOLUMES_H

#include <string>
#include <vector>

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
class LayerVolumes
{
public:
	LayerVolumes();
	~LayerVolumes() {}

	bool createGeoVolumes(const MeshLib::Mesh &mesh, const std::vector<std::string> &raster_paths, double noDataReplacementValue = 0.0);

	MeshLib::Mesh* getMesh() const { return _mesh; }

	bool addGeometry(GeoLib::GEOObjects &geo_objects) const;

	std::vector<MeshLib::Node> getAttributePoints() { return _attribute_points; }

private:
	void addLayerToMesh(const MeshLib::Mesh &mesh_layer, unsigned layer_id);

	void addLayerBoundaries(const MeshLib::Mesh &layer, std::size_t nLayers);

	void removeCongruentElements(std::size_t nLayers, std::size_t nElementsPerLayer);

	bool allRastersExist(const std::vector<std::string> &raster_paths) const;

	void cleanUpOnError();

	const double _invalid_value;
	std::vector<MeshLib::Node*> _nodes;
	std::vector<MeshLib::Element*> _elements;
	std::vector<MeshLib::Node> _attribute_points;
	MeshLib::Mesh* _mesh;
};

#endif //LAYERVOLUMES_H
