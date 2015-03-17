/**
 * \file
 * \author Karsten Rink
 * \date   2013-04-04
 * \brief  Definition of the ElementValueModification class
 *
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MESHVALUEMODIFICATION_H
#define MESHVALUEMODIFICATION_H

#include <vector>

#include <boost/optional.hpp>

#include "MeshLib/MeshEnums.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/PropertyVector.h"

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

	/// Sets new value for all elements having the given element type
	/// Returns the number of elements having the given element type
	static unsigned setByElementType(MeshLib::Mesh &mesh, MeshElemType ele_type, unsigned new_value);

private:
	/// Returns sorted values of properties within the PropertyVector
	/// These values are stored in a vector.
	template <typename T>
	static std::vector<T> getSortedPropertyValues(
		MeshLib::PropertyVector<T> const& property_vector)
	{
		std::vector<T> value_mapping;
		const std::size_t n_property_values(property_vector.size());
		for (std::size_t i=0; i<n_property_values; ++i) {
			bool exists(false);
			T const& value (property_vector[i]);
			std::size_t const size(value_mapping.size());
			for (unsigned j=0; j<size; ++j) {
				if (value == value_mapping[j]) {
					exists = true;
					break;
				}
			}
			if (!exists)
				value_mapping.push_back(value);
		}

		std::sort(value_mapping.begin(), value_mapping.end());
		return value_mapping;
	}

	/// Returns the values of elements within the mesh
	static std::vector<unsigned> getMeshValues(const MeshLib::Mesh &mesh);
};

} // end namespace MeshLib

#endif //MESHVALUEMODIFICATION_H
