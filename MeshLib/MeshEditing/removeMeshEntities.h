/**
 * \file
 * \author Karsten Rink
 * \date   2013-04-04
 * \brief  Definition of the removeMeshEntities
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef REMOVEMESHENTITIES_H
#define REMOVEMESHENTITIES_H

#include <vector>
#include "MeshEnums.h"

namespace MeshLib {

// forward declarations
class Mesh;

	/// Removes the mesh nodes (and connected elements) given in the nodes-list from the mesh.
	MeshLib::Mesh* removeMeshNodes(MeshLib::Mesh const*const mesh, const std::vector<std::size_t> &nodes);

	/// Removes elements of the given type t from a mesh
	MeshLib::Mesh* removeMeshElements(const MeshLib::Mesh &mesh, MeshElemType t);

} // end namespace MeshLib

#endif //REMOVEMESHENTITIES_H
