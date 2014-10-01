/**
 * \file   MeshSurfaceExtraction.h
 * \author Karsten Rink
 * \date   2013-04-04
 * \brief  Definition of the MeshSurfaceExtraction class
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MESHSURFACEEXTRACTION_H
#define MESHSURFACEEXTRACTION_H

#include <cstddef>
#include <vector>

#include "MathLib/Vector3.h"

namespace GeoLib {
	class PointWithID;
}

namespace MeshLib {
// forward declarations
class Mesh;
class Element;
class Node;

/**
 * \brief A set of tools concerned with extracting nodes and elements from a mesh surface
 */
class MeshSurfaceExtraction
{
public:
	/// Returns a vector of the areas assigned to each node on a surface mesh.
	static std::vector<double> getSurfaceAreaForNodes(const MeshLib::Mesh &mesh);

	/// Returns the surface nodes of a layered mesh.
	static std::vector<GeoLib::PointWithID*> getSurfaceNodes(const MeshLib::Mesh &mesh, const MathLib::Vector3 &dir, double angle);

	/**
	 * Returns the 2d-element mesh representing the surface of the given layered mesh.
	 * \param mesh                The original mesh
	 * \param dir                 The direction in which face normals have to point to be considered surface elements
	 * \param angle               The angle of the allowed deviation from the given direction (0 <= angle <= 90 degrees)
	 * \param keepOriginalNodeIds If true, ids of mesh nodes are set to ids in original mesh, otherwise node ids are reset (as usual when creating a mesh)
	 * \return                    A 2D mesh representing the surface in direction dir
	 */
	static MeshLib::Mesh* getMeshSurface(const MeshLib::Mesh &mesh, const MathLib::Vector3 &dir, double angle, bool keepOriginalNodeIds = false);

private:
	/// Functionality needed for getSurfaceNodes() and getMeshSurface()
	static void get2DSurfaceElements(const std::vector<MeshLib::Element*> &all_elements, std::vector<MeshLib::Element*> &sfc_elements, const MathLib::Vector3 &dir, double angle, unsigned mesh_dimension);

	/// Functionality needed for getSurfaceNodes() and getMeshSurface()
	static void get2DSurfaceNodes(std::vector<MeshLib::Node*> &sfc_nodes, std::size_t n_all_nodes, const std::vector<MeshLib::Element*> &sfc_elements, std::vector<std::size_t> &node_id_map);
};

} // end namespace MeshLib

#endif //MESHSURFACEEXTRACTION_H
