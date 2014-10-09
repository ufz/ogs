/**
 * \file
 * \author Karsten Rink
 * \date   2013-10-28
 * \brief  Definition of the MeshInformation class
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MESHINFORMATION_H
#define MESHINFORMATION_H

#include <array>
#include "GeoLib/AABB.h"
#include "MeshLib/Node.h"

namespace MeshLib
{
class Mesh;

/**
 * \brief A set of tools for extracting information from a mesh
 */
class MeshInformation
{
public:
	/// Returns the smallest and largest MaterialID of the mesh.
	static const std::pair<unsigned, unsigned> getValueBounds(const MeshLib::Mesh &mesh);

	/// Returns the bounding box of the mesh.
	static const GeoLib::AABB<MeshLib::Node> getBoundingBox(const MeshLib::Mesh &mesh);

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
	static const std::array<unsigned, 7> getNumberOfElementTypes(const MeshLib::Mesh &mesh);


};

}

#endif //MESHINFORMATION_H
