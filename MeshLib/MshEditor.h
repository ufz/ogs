/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
 *
 * \file MshEditor.h
 *
 * Created on 2011-06-15 by Karsten Rink
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

/**
 * \brief A set of tools for manipulating existing meshes
 */
class MshEditor
{
public:
	MshEditor() {}
	~MshEditor() {}

	/// Returns the area assigned to each node on a surface mesh.
	static void getSurfaceAreaForNodes(const MeshLib::Mesh* mesh, std::vector<double> &node_area_vec);

	/// Removes the mesh nodes (and connected elements) given in the nodes-list from the mesh.
	static MeshLib::Mesh* removeMeshNodes(MeshLib::Mesh* mesh, const std::vector<size_t> &nodes);

	/// Returns the surface nodes of a layered mesh.
	static std::vector<GeoLib::PointWithID*> getSurfaceNodes(const MeshLib::Mesh &mesh);

	/// Returns the 2d-element mesh representing the surface of the given layered mesh.
	static MeshLib::Mesh* getMeshSurface(const MeshLib::Mesh &mesh, const double* dir);

private:
};

} // end namespace MeshLib

#endif //MSHEDITOR_H
