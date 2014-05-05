/**
 * \file
 * \author Karsten Rink
 * \date   2013-04-04
 * \brief  Definition of the MeshSurfaceExtraction class
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MESHSURFACEEXTRACTION_H
#define MESHSURFACEEXTRACTION_H

#include <cstddef>
#include <vector>

#include "Vector3.h"

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
	/// Returns the area assigned to each node on a surface mesh.
	static void getSurfaceAreaForNodes(const MeshLib::Mesh* mesh, std::vector<double> &node_area_vec);

	/// Returns the surface nodes of a layered mesh.
	static std::vector<GeoLib::PointWithID*> getSurfaceNodes(const MeshLib::Mesh &mesh, const MathLib::Vector3 &dir);

	/// Returns the 2d-element mesh representing the surface of the given layered mesh.
	static MeshLib::Mesh* getMeshSurface(const MeshLib::Mesh &mesh, const MathLib::Vector3 &dir);

private:
	/// Functionality needed for getSurfaceNodes() and getMeshSurface()
	static void get2DSurfaceElements(const std::vector<MeshLib::Element*> &all_elements, std::vector<MeshLib::Element*> &sfc_elements, const MathLib::Vector3 &dir, unsigned mesh_dimension);

	/// Functionality needed for getSurfaceNodes() and getMeshSurface()
	static void get2DSurfaceNodes(const std::vector<MeshLib::Node*> &all_nodes, std::vector<MeshLib::Node*> &sfc_nodes, const std::vector<MeshLib::Element*> &sfc_elements, std::vector<unsigned> &node_id_map);
};

} // end namespace MeshLib

#endif //MESHSURFACEEXTRACTION_H
