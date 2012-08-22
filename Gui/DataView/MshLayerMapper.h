/**
 * \file MshLayerMapper.h
 * 01/11/2010 KR Initial implementation
 */

#ifndef MSHLAYERMAPPER_H
#define MSHLAYERMAPPER_H

#include <string>

class QImage;

namespace MeshLib
{
class CFEMesh;
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
	static MeshLib::CFEMesh* CreateLayers(const MeshLib::CFEMesh* mesh,
	                                      size_t nLayers,
	                                      double thickness);

	/// Maps the z-values of nodes in the designated layer of the given mesh according to the given raster.
	static int LayerMapping(MeshLib::CFEMesh* msh,
	                                      const std::string &rasterfile,
	                                      const size_t nLayers,
	                                      const size_t layer_id,
	                                      bool removeNoDataValues = false);

	/// Blends a mesh with the surface given by dem_raster. Nodes and elements above the surface are either removed or adapted to fit the surface.
	static MeshLib::CFEMesh* blendLayersWithSurface(MeshLib::CFEMesh* mesh, const size_t nLayers, const std::string &dem_raster);

private:
	/// Checks if the given mesh is within the dimensions given by xDim and yDim.
	static bool meshFitsImage(const MeshLib::CFEMesh* msh,
	                          const std::pair<double, double> &xDim,
	                          const std::pair<double, double> &yDim);
};

#endif //MSHLAYERMAPPER_H
