/**
 * \file
 * \author Karsten Rink
 * \date   2013-04-04
 * \brief  Definition of the ElementValueModification class
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MESHVALUEMODIFICATION_H
#define MESHVALUEMODIFICATION_H

#include <vector>

namespace MeshLib {
// forward declarations
class Mesh;

/**
 * \brief A set of methods for manipulating mesh element values
 */
class ElementValueModification
{
public:
	/// Reduces the values assigned the elements of mesh to the smallest possible range.
	/// Returns the number of different values.
	static unsigned condense(MeshLib::Mesh &mesh);

	/// Replaces for all elements of mesh with the value old_value with new_value if possible.
	/// Returns true if successful or false if the value is already taken.
	static bool replace(MeshLib::Mesh &mesh, unsigned old_value, unsigned new_value, bool replace_if_exists = false);

private:
	/// Returns the values of elements within the mesh
	static std::vector<unsigned> getMeshValues(const MeshLib::Mesh &mesh);
};

} // end namespace MeshLib

#endif //MESHVALUEMODIFICATION_H
