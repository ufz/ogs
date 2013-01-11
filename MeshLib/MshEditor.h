/**
 * \file
 * \author Karsten Rink
 * \date   2011-06-15
 * \brief  Definition of the MshEditor class
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MSHEDITOR_H
#define MSHEDITOR_H

#include <cstddef>
#include <vector>

#include "Point.h"

namespace GeoLib {
	class PointWithID;
}

namespace MeshLib {
// forward declarations
class Mesh;
class Element;
class Node;

/**
 * \brief A set of tools for manipulating existing meshes
 */
class MshEditor
{
public:
	/// Returns the area assigned to each node on a surface mesh.
	static void getSurfaceAreaForNodes(const MeshLib::Mesh* mesh, std::vector<double> &node_area_vec);

	/// Removes the mesh nodes (and connected elements) given in the nodes-list from the mesh.
	static MeshLib::Mesh* removeMeshNodes(MeshLib::Mesh* mesh, const std::vector<std::size_t> &nodes);

	/// Returns the surface nodes of a layered mesh.
	static std::vector<GeoLib::PointWithID*> getSurfaceNodes(const MeshLib::Mesh &mesh, const double* dir = NULL);

	/// Returns the 2d-element mesh representing the surface of the given layered mesh.
	static MeshLib::Mesh* getMeshSurface(const MeshLib::Mesh &mesh, const double* dir = NULL);

private:
	/// Functionality needed for getSurfaceNodes() and getMeshSurface()
	static void get2DSurfaceElements(const std::vector<MeshLib::Element*> &all_elements, std::vector<MeshLib::Element*> &sfc_elements, const double* dir, unsigned mesh_dimension);

	/// Functionality needed for getSurfaceNodes() and getMeshSurface()
	static void get2DSurfaceNodes(const std::vector<MeshLib::Node*> &all_nodes, std::vector<MeshLib::Node*> &sfc_nodes, const std::vector<MeshLib::Element*> &sfc_elements, std::vector<unsigned> &node_id_map);
};

} // end namespace MeshLib

#endif //MSHEDITOR_H
