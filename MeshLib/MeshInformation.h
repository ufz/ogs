/**
 * \file
 * \author Karsten Rink
 * \date   2013-10-28
 * \brief  Definition of the MeshInformation class
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MESHINFORMATION_H
#define MESHINFORMATION_H

#include <array>

#include "AABB.h"
#include "Node.h"

namespace MeshLib {
	class Mesh;
}

/**
 * \brief A set of tools for extracting information from a mesh
 */
class MeshInformation
{
public:
	/// Returns the smallest and largest MaterialID of the mesh.
	static const std::pair<unsigned, unsigned> getValueBounds(MeshLib::Mesh const*const mesh);

	/// Returns the bounding box of the mesh.
	static const GeoLib::AABB<MeshLib::Node> getBoundingBox(MeshLib::Mesh const*const mesh);

	/** 
	 * Returns an array with the number of elements of each type in the given mesh.
	 * On completion, n_element_types array contains the number of elements of each of the seven 
	 * supported types. The index to element type conversion is this:
	 *		0: #lines
	 *		1: #triangles
	 *		2: #quads
	 *		3: #tetrahedra
	 *		4: #hexahedra
	 *		5: #pyramids
	 *		6: #prisms
	 */
	static void getNumberOfElementTypes(MeshLib::Mesh const*const mesh, std::array<unsigned, 7> &n_element_types);


};

#endif //MESHINFORMATION_H
