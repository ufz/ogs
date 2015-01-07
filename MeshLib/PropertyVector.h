/**
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROPERTYVECTOR_H_
#define PROPERTYVECTOR_H_

#include <algorithm>
#include <vector>
#include <type_traits>
#include <memory>

namespace MeshLib
{

/// Class template PropertyVector is a std::vector with template parameter
/// PROP_VAL_TYPE. The reason for the derivation of std::vector is
/// the template specialisation for pointer types below.
/// \tparam PROP_VAL_TYPE typical this is a scalar, a vector or a matrix
template <typename PROP_VAL_TYPE>
class PropertyVector : public std::vector<PROP_VAL_TYPE>
{
friend class Properties;
protected:
	PropertyVector()
		: std::vector<PROP_VAL_TYPE>()
	{}

	explicit PropertyVector(std::size_t size)
		: std::vector<PROP_VAL_TYPE>(size)
	{}
};

/// Class template PropertyVector is a std::vector with template parameter
/// T, where T is a pointer type.
/// The behaviour has changed for the constructor, destructor and the operator[].
/// The user has to provide the size and an item to group mapping for construction.
/// The destructor takes care to delete the entries of the vector.
/// The operator[] uses an item-to-group property map to access the
/// correct property.
/// \tparam T pointer type, the type the type points to is typical a scalar,
/// a vector or a matrix type
template <typename T>
class PropertyVector<T*> : public std::vector<T*>
{
friend class Properties;
protected:
	/// @param n_prop_groups number of different property values
	/// @param item2group_mapping Class Mesh has a mapping from the mesh items
	/// (Node or Element) to an index (position in the data structure).
	/// The vector item2group_mapping must have the same number of entries as
	/// the above mapping and the values have to be in the range
	/// \f$[0, \text{n_prop_groups})\f$.
	PropertyVector(std::size_t n_prop_groups,
		std::vector<std::size_t> const& item2group_mapping)
		: std::vector<T*>(n_prop_groups),
		_item2group_mapping(item2group_mapping)
	{}

public:
	/// Destructor ensures the deletion of the heap-constructed objects.
	~PropertyVector()
	{
		std::for_each(
			this->begin(), this->end(), std::default_delete<T>()
		);
	}

	/// The operator[] uses the item to group property map to access to the
	/// correct property value/object.
	T* const& operator[](std::size_t id) const
	{
		return (*static_cast<std::vector<T*> const*>(this))[_item2group_mapping[id]];
	}

private:
	std::vector<std::size_t> _item2group_mapping;
};

} // end namespace MeshLib

#endif
