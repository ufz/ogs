/**
 * \file
 * \brief  Implementation of the class Properties that implements a container of
 *         properties.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "Properties.h"

namespace MeshLib
{

void Properties::removePropertyVector(std::string const& name)
{
    std::map<std::string, PropertyVectorBase*>::const_iterator it(
        _properties.find(name)
    );
    if (it == _properties.end()) {
        WARN("A property of the name \"%s\" does not exist.",
            name.c_str());
        return;
    }
    delete it->second;
    _properties.erase(it);
}

bool Properties::hasPropertyVector(std::string const& name) const
{
    std::map<std::string, PropertyVectorBase*>::const_iterator it(
        _properties.find(name)
    );
    if (it == _properties.end()) {
        return false;
    }
    return true;
}

std::vector<std::string> Properties::getPropertyVectorNames() const
{
    std::vector<std::string> names;
    for (auto p : _properties)
        names.push_back(p.first);
    return names;
}

Properties Properties::excludeCopyProperties(
    std::vector<std::size_t> const& exclude_elem_ids,
    std::vector<std::size_t> const& exclude_node_ids) const
{
    Properties exclude_copy;
    for (auto property_vector : _properties) {
        if (property_vector.second->getMeshItemType() == MeshItemType::Cell) {
            exclude_copy._properties.insert(
                std::make_pair(property_vector.first,
                property_vector.second->clone(exclude_elem_ids))
            );
        }
        else if (property_vector.second->getMeshItemType() == MeshItemType::Node)
        {
            exclude_copy._properties.insert(
                std::make_pair(property_vector.first,
                property_vector.second->clone(exclude_node_ids))
            );
        }
    }
    return exclude_copy;
}

Properties Properties::excludeCopyProperties(
    std::vector<MeshItemType> const& exclude_mesh_item_types) const
{
    Properties new_properties;
    for (auto property_vector : _properties) {
        if (std::find(exclude_mesh_item_types.begin(),
                      exclude_mesh_item_types.end(),
                      property_vector.second->getMeshItemType()) !=
            exclude_mesh_item_types.end())
            continue;

        std::vector<std::size_t> const exclude_positions{};
        new_properties._properties.insert(
            std::make_pair(property_vector.first,
            property_vector.second->clone(exclude_positions))
        );
    }
    return new_properties;
}

Properties::Properties(Properties const& properties)
    : _properties(properties._properties)
{
    std::vector<std::size_t> exclude_positions;
    for (auto it(_properties.begin()); it != _properties.end(); ++it) {
        PropertyVectorBase *t(it->second->clone(exclude_positions));
        it->second = t;
    }
}

Properties::~Properties()
{
    for (auto property_vector : _properties) {
        delete property_vector.second;
    }
}

} // end namespace MeshLib

