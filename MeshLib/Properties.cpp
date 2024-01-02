/**
 * \file
 * \brief  Implementation of the class Properties that implements a container of
 *         properties.
 *
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "Properties.h"

namespace MeshLib
{
void Properties::removePropertyVector(std::string_view name)
{
    auto it(_properties.find(std::string(name)));
    if (it == _properties.end())
    {
        WARN("A property of the name '{:s}' does not exist.", name);
        return;
    }
    delete it->second;
    _properties.erase(it);
}

bool Properties::hasPropertyVector(std::string_view name) const
{
    return _properties.find(std::string(name)) != _properties.end();
}

std::vector<std::string> Properties::getPropertyVectorNames() const
{
    std::vector<std::string> names;
    std::transform(_properties.begin(), _properties.end(),
                   std::back_inserter(names),
                   [](auto const& pair) { return std::string(pair.first); });
    return names;
}

std::vector<std::string> Properties::getPropertyVectorNames(
    MeshLib::MeshItemType t) const
{
    std::vector<std::string> names;
    for (auto p : _properties)
    {
        if (p.second->getMeshItemType() == t)
        {
            names.push_back(std::string(p.first));
        }
    }
    return names;
}

Properties Properties::excludeCopyProperties(
    std::vector<std::size_t> const& exclude_elem_ids,
    std::vector<std::size_t> const& exclude_node_ids) const
{
    Properties exclude_copy;
    for (auto name_vector_pair : _properties)
    {
        if (name_vector_pair.second->getMeshItemType() == MeshItemType::Cell)
        {
            exclude_copy._properties.insert(std::make_pair(
                name_vector_pair.first,
                name_vector_pair.second->clone(exclude_elem_ids)));
        }
        else if (name_vector_pair.second->getMeshItemType() ==
                 MeshItemType::Node)
        {
            exclude_copy._properties.insert(std::make_pair(
                name_vector_pair.first,
                name_vector_pair.second->clone(exclude_node_ids)));
        }
    }
    return exclude_copy;
}

Properties Properties::excludeCopyProperties(
    std::vector<MeshItemType> const& exclude_mesh_item_types) const
{
    Properties new_properties;
    for (auto name_vector_pair : _properties)
    {
        if (std::find(exclude_mesh_item_types.begin(),
                      exclude_mesh_item_types.end(),
                      name_vector_pair.second->getMeshItemType()) !=
            exclude_mesh_item_types.end())
        {
            continue;
        }

        std::vector<std::size_t> const exclude_positions{};
        new_properties._properties.insert(
            std::make_pair(name_vector_pair.first,
                           name_vector_pair.second->clone(exclude_positions)));
    }
    return new_properties;
}

Properties::Properties(Properties const& properties)
    : _properties(properties._properties)
{
    for (auto& name_vector_pair : _properties)
    {
        PropertyVectorBase* t(name_vector_pair.second->clone({}));
        name_vector_pair.second = t;
    }
}

Properties& Properties::operator=(Properties const& properties)
{
    if (&properties == this)
    {
        return *this;
    }

    _properties = properties._properties;
    std::vector<std::size_t> exclude_positions;
    for (auto& name_vector_pair : _properties)
    {
        PropertyVectorBase* t(
            name_vector_pair.second->clone(exclude_positions));
        name_vector_pair.second = t;
    }

    return *this;
}

Properties::~Properties()
{
    for (auto name_vector_pair : _properties)
    {
        delete name_vector_pair.second;
    }
}

std::map<std::string, PropertyVectorBase*>::const_iterator Properties::begin()
    const
{
    return _properties.cbegin();
}

std::map<std::string, PropertyVectorBase*>::const_iterator Properties::end()
    const
{
    return _properties.cend();
}

std::map<std::string, PropertyVectorBase*>::iterator Properties::begin()
{
    return _properties.begin();
}

std::map<std::string, PropertyVectorBase*>::iterator Properties::end()
{
    return _properties.end();
}

std::map<std::string, PropertyVectorBase*>::size_type Properties::size() const
{
    return _properties.size();
}

std::map<std::string, PropertyVectorBase*>::size_type Properties::size(
    MeshItemType const mesh_item_type) const
{
    return count_if(begin(), end(),
                    [&](auto const p)
                    { return p.second->getMeshItemType() == mesh_item_type; });
}

}  // end namespace MeshLib
