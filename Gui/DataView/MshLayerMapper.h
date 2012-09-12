/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
 *
 * \file MshLayerMapper.h
 *
 * Created on 2010-11-01 by Karsten Rink
 */

#ifndef MSHLAYERMAPPER_H
#define MSHLAYERMAPPER_H

#include <string>

class QImage;

namespace MeshLib {
	class Mesh;
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
	static MeshLib::Mesh* CreateLayers(const MeshLib::Mesh* mesh, std::size_t nLayers, double thickness);

	/// Maps the z-values of nodes in the designated layer of the given mesh according to the given raster.
	static int LayerMapping(MeshLib::Mesh* msh, const std::string &rasterfile,
                            const std::size_t nLayers, const std::size_t layer_id, bool removeNoDataValues = false);

	/// Blends a mesh with the surface given by dem_raster. Nodes and elements above the surface are either removed or adapted to fit the surface.
	static MeshLib::Mesh* blendLayersWithSurface(MeshLib::Mesh* mesh, const std::size_t nLayers, const std::string &dem_raster);

private:
	/// Checks if the given mesh is within the dimensions given by xDim and yDim.
	static bool meshFitsImage(const MeshLib::Mesh* msh,
	                          const std::pair<double, double> &xDim,
	                          const std::pair<double, double> &yDim);
};

#endif //MSHLAYERMAPPER_H
