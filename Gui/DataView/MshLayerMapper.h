/**
 * \file
 * \author Karsten Rink
 * \date   2010-11-01
 * \brief  Definition of the MshLayerMapper class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MSHLAYERMAPPER_H
#define MSHLAYERMAPPER_H

#include <string>

class QImage;

namespace MeshLib {
	class Mesh;
	class Node;
}

/**
 * \brief Manipulating and adding layers to an existing mesh
 */
class MshLayerMapper
{
public:
	MshLayerMapper() {}
	~MshLayerMapper() {}

	/**
	 * Based on a triangle-or quad mesh this method creates a 3D mesh with with a given number of prism- or hex-layers
	 * \param mesh The triangle/quad mesh that is the basis for the new prism/hex mesh
	 * \param nLayers The number of layers of prism/hex elements that will be extruded from the triangle/quad elements of the original mesh
	 * \param thickness The thickness of each of these newly added layers
	 * \return A mesh with the requested number of layers of prism/hex elements
	 */
	static MeshLib::Mesh* CreateLayers(const MeshLib::Mesh &mesh, const std::vector<float> &layer_thickness_vector);

	/**
	 * Maps the z-values of nodes in the designated layer of the given mesh according to the given raster.
	 * Note: This only results in a valid mesh if the layers don't intersect each other.
	 */
	static int LayerMapping(MeshLib::Mesh &mesh, const std::string &rasterfile,
                            const unsigned nLayers, const unsigned layer_id, double noDataReplacementValue);

	/**
	 * Blends a mesh with the surface given by dem_raster. Nodes and elements above the surface are either removed or adapted to fit the surface.
	 * Note: It is unlikely but possible that the new nodes vector contains (very few) nodes that are not part of any element. This problem is
	 * remedied at the end of method upon creating the actual mesh from the new node- and element-vector as the mesh-constructor checks for such
	 * nodes and removes them. This note is just to call this issue to attention in case this methods is changed.
	 */
	static MeshLib::Mesh* blendLayersWithSurface(MeshLib::Mesh* mesh, const unsigned nLayers, const std::string &dem_raster);

private:
	/// Checks if the given mesh is within the dimensions given by xDim and yDim.
	static bool isNodeOnRaster(const MeshLib::Node &node,
	                           const std::pair<double, double> &xDim,
	                           const std::pair<double, double> &yDim);
};

#endif //MSHLAYERMAPPER_H
