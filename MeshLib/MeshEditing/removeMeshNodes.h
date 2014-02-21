/**
 * \file
 * \author Karsten Rink
 * \date   2013-04-04
 * \brief  Definition of the removeMeshNodes
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef REMOVEMESHNODES_H
#define REMOVEMESHNODES_H

#include <vector>

namespace MeshLib {

// forward declarations
class Mesh;

	/// Removes the mesh nodes (and connected elements) given in the nodes-list from the mesh.
	/// Warning: this function actually modifies the mesh, it might make sense to copy the mesh before using this function.
	void removeMeshNodes(MeshLib::Mesh &mesh, const std::vector<std::size_t> &nodes);

} // end namespace MeshLib

#endif //REMOVEMESHNODES_H
