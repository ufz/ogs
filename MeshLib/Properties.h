/**
 * \file
 * \brief  Definition of the class Properties that implements a container of
 *         properties.
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROPERTIES_H_
#define PROPERTIES_H_

#include <cstdlib>
#include <string>
#include <map>

#include <boost/any.hpp>
#include <boost/optional.hpp>

#include "logog/include/logog.hpp"

#include "Location.h"

#include "PropertyVector.h"

namespace MeshLib
{

/// Properties associated to mesh items (nodes or elements).
class Properties
{
public:
	template <typename T>
	boost::optional<PropertyVector<T> &>
	newProperty(std::string const& name, MeshItemType mesh_item_type)
	{
		PropertyKeyType property_key(name, mesh_item_type);
		std::map<PropertyKeyType, boost::any>::const_iterator it(
			_properties.find(property_key)
		);
		if (it != _properties.end()) {
			WARN("A property of the name \"%s\" already assigned to the mesh.",
				name.c_str());
			return boost::optional<PropertyVector<T> &>();
		}
		auto entry_info(
			_properties.insert(
				std::pair<PropertyKeyType, boost::any>(
					property_key, boost::any(PropertyVector<T>())
				)
			)
		);
		return boost::optional<PropertyVector<T> &>(
			boost::any_cast<PropertyVector<T> &>((entry_info.first)->second)
			);
	}

	template <typename T>
	boost::optional<PropertyVector<T> &>
	newProperty(std::string const& name,
		std::size_t n_prop_groups,
		std::vector<std::size_t> const& item2group_mapping,
		MeshItemType mesh_item_type)
	{
		PropertyKeyType property_key(name, mesh_item_type);
		std::map<PropertyKeyType, boost::any>::const_iterator it(
			_properties.find(property_key)
		);
		if (it != _properties.end()) {
			WARN("A property of the name \"%s\" already assigned to the mesh.",
				name.c_str());
			return boost::optional<PropertyVector<T> &>();
		}
		auto entry_info(
			_properties.insert(
				std::pair<PropertyKeyType, boost::any>(
					property_key,
					boost::any(PropertyVector<T>(n_prop_groups, item2group_mapping))
				)
			)
		);
		return boost::optional<PropertyVector<T> &>(
			boost::any_cast<PropertyVector<T> &>(
				(entry_info.first)->second)
			);
	}

	/// Method to get a vector of property values.
	template <typename T>
	boost::optional<PropertyVector<T> const&>
	getProperty(std::string const& name,
		MeshItemType mesh_item_type)
	{
		PropertyKeyType property_key(name, mesh_item_type);
		std::map<PropertyKeyType, boost::any>::const_iterator it(
			_properties.find(property_key)
		);
		if (it != _properties.end()) {
			try {
				return boost::optional<PropertyVector<T> const&>(
						boost::any_cast<PropertyVector<T> const&>(it->second)
					);
			} catch (boost::bad_any_cast const&) {
				ERR("A property with the desired data type is not available.");
				return boost::optional<PropertyVector<T> const&>();
			}
		} else {
			return boost::optional<PropertyVector<T> const&>();
		}
	}

	/// Method to store a vector of property values assigned to a property name.
	/// Since the implementation makes no assumption about the number of data
	/// items stored within the vector, it is possible either to use a small
	/// number of properties where each particular property can be assigned to
	/// several mesh items. In contrast to this it is possible to have a
	/// separate value for each mesh item.
	/// The user has to ensure the correct usage of the vector later on.
	template <typename T>
	void addProperty(std::string const& name, PropertyVector<T> * property,
		MeshItemType mesh_item_type)
	{
		PropertyKeyType property_key(name, mesh_item_type);
		std::map<PropertyKeyType, boost::any>::const_iterator it(
			_properties.find(property_key)
		);
		if (it != _properties.end()) {
			WARN("A property of the name \"%s\" already assigned to the mesh.",
				name.c_str());
			return;
		}
		_properties[property_key] = boost::any(property);
	}

	void removeProperty(std::string const& name,
		MeshItemType mesh_item_type)
	{
		PropertyKeyType property_key(name, mesh_item_type);
		std::map<PropertyKeyType, boost::any>::const_iterator it(
			_properties.find(property_key)
		);
		if (it == _properties.end()) {
			WARN("A property of the name \"%s\" does not exist.",
				name.c_str());
			return;
		}
		_properties.erase(it);
	}

private:
	struct PropertyKeyType
	{
		PropertyKeyType(std::string const& n, MeshItemType t)
			: name(n), mesh_item_type(t)
		{}

		std::string name;
		MeshItemType mesh_item_type;

		bool operator<(PropertyKeyType const& other) const
		{
			if (name.compare(other.name) == 0) {
				return mesh_item_type < other.mesh_item_type;
			}
			return name.compare(other.name) < 0;
		}
	};

	/// A mapping from property's name to the stored object of any type.
	/// See addProperty() and getProperty() documentation.
	std::map<PropertyKeyType, boost::any> _properties;
}; // end class

} // end namespace MeshLib

#endif

