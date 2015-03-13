/**
 *
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
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

#include "BaseLib/excludeObjectCopy.h"

namespace MeshLib
{

template <typename T>
class PropertyVector;

class PropertyVectorBase
{
public:
	virtual PropertyVectorBase* clone(
		std::vector<std::size_t> const& exclude_positions
	) const = 0;
	virtual ~PropertyVectorBase() = default;
};

/// Class template PropertyVector is a std::vector with template parameter
/// PROP_VAL_TYPE. The reason for the derivation of std::vector is
/// the template specialisation for pointer types below.
/// \tparam PROP_VAL_TYPE typical this is a scalar, a vector or a matrix
template <typename PROP_VAL_TYPE>
class PropertyVector : public std::vector<PROP_VAL_TYPE>,
	public PropertyVectorBase
{
friend class Properties;

public:
	std::size_t getTupleSize() const { return _tuple_size; }
	MeshItemType getMeshItemType() const { return _mesh_item_type; }
	std::string const& getPropertyName() const { return _property_name; }

	PropertyVectorBase* clone(std::vector<std::size_t> const& exclude_positions) const
	{
		PropertyVector<PROP_VAL_TYPE> *t(new PropertyVector<PROP_VAL_TYPE>(_property_name,
			_mesh_item_type, _tuple_size));
		BaseLib::excludeObjectCopy(*this, exclude_positions, *t);
		return t;
	}

protected:
	/// @brief The constructor taking meta information for the data.
	/// @param property_name a string describing the property
	/// @param mesh_item_type the values of the property are either assigned to
	/// nodes or cells (see enumeration @MeshItemType) (default:
	/// MeshItemType::Cell)
	/// @param tuple_size the number of elements of a tuple (default: 1)
	explicit PropertyVector(std::string const& property_name,
		MeshItemType mesh_item_type = MeshItemType::Cell,
		std::size_t tuple_size = 1)
		: std::vector<PROP_VAL_TYPE>(), _tuple_size(tuple_size),
		_mesh_item_type(mesh_item_type), _property_name(property_name)
	{}

	/// @brief The constructor taking meta information for the data.
	/// @param n_property_values number of property values (value can be a tuple
	/// with several entries)
	/// @param property_name a string describing the property
	/// @param mesh_item_type the values of the property are either assigned to
	/// nodes or cells (see enumeration @MeshItemType) (default:
	/// MeshItemType::Cell)
	/// @param tuple_size the number of elements of a tuple (default: 1)
	PropertyVector(std::size_t n_property_values,
		std::string const& property_name,
		MeshItemType mesh_item_type = MeshItemType::Cell,
		std::size_t tuple_size = 1)
		: std::vector<PROP_VAL_TYPE>(n_property_values*tuple_size)
	{}

	std::size_t const _tuple_size;
	MeshItemType const _mesh_item_type;
	std::string const _property_name;
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
class PropertyVector<T*> : public std::vector<T*>,
	public PropertyVectorBase
{
friend class Properties;
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

	T* & operator[](std::size_t id)
	{
		return (*static_cast<std::vector<T*>*>(this))[_item2group_mapping[id]];
	}

	std::size_t getTupleSize() const { return _tuple_size; }
	MeshItemType getMeshItemType() const { return _mesh_item_type; }
	std::string const& getPropertyName() const { return _property_name; }

	PropertyVectorBase* clone(std::vector<std::size_t> const& exclude_positions) const
	{
		std::vector<std::size_t> item2group_mapping(
			BaseLib::excludeObjectCopy(_item2group_mapping, exclude_positions)
		);
		PropertyVector<T*> *t (new PropertyVector<T*>(this->size()/_tuple_size,
			item2group_mapping, _property_name, _mesh_item_type, _tuple_size));
		for (std::size_t j(0); j<item2group_mapping.size(); j++) {
			t->operator[](j) = (*static_cast<std::vector<T*> const*>(this))[item2group_mapping[j]];
		}
		return t;
	}

protected:
	/// @brief The constructor taking meta information for the data.
	/// @param n_prop_groups number of different property values
	/// @param item2group_mapping Class Mesh has a mapping from the mesh items
	/// (Node or Element) to an index (position in the data structure).
	/// The vector item2group_mapping must have the same number of entries as
	/// the above mapping and the values have to be in the range
	/// \f$[0, \text{n_prop_groups})\f$.
	/// @param property_name a string describing the property
	/// @param mesh_item_type the values of the property are either assigned to
	/// nodes or cells (see enumeration @MeshItemType) (default:
	/// MeshItemType::Cell)
	/// @param tuple_size the number of elements of a tuple (default: 1)
	PropertyVector(std::size_t n_prop_groups,
		std::vector<std::size_t> const& item2group_mapping,
		std::string const& property_name,
		MeshItemType mesh_item_type = MeshItemType::Cell,
		std::size_t tuple_size = 1)
		: std::vector<T*>(n_prop_groups * tuple_size),
		_tuple_size(tuple_size),
		_mesh_item_type(mesh_item_type),
		_property_name(property_name),
		_item2group_mapping(item2group_mapping)
	{}

protected:
	std::size_t const _tuple_size;
	MeshItemType const _mesh_item_type;
	std::string const _property_name;

private:
	T* at(std::size_t);
	std::vector<std::size_t> _item2group_mapping;
};

} // end namespace MeshLib

#endif
