/**
 * \file
 * \brief  Implementation of the class Properties that implements a container of
 *         properties.
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "Properties.h"

namespace MeshLib
{

void Properties::removePropertyVector(std::string const& name)
{
    std::map<std::string, PropertyVectorBase*>::const_iterator it(
        properties_.find(name)
    );
    if (it == properties_.end()) {
        WARN("A property of the name '{:s}' does not exist.", name);
        return;
    }
    delete it->second;
    properties_.erase(it);
}

bool Properties::hasPropertyVector(std::string const& name) const
{
    return properties_.find(name) != properties_.end();
}

std::vector<std::string> Properties::getPropertyVectorNames() const
{
    std::vector<std::string> names;
    std::transform(properties_.begin(), properties_.end(),
                   std::back_inserter(names),
                   [](auto const& pair) { return pair.first; });
    return names;
}

std::vector<std::string> Properties::getPropertyVectorNames(
    MeshLib::MeshItemType t) const
{
    std::vector<std::string> names;
    for (auto p : properties_)
    {
        if (p.second->getMeshItemType() == t)
        {
            names.push_back(p.first);
        }
    }
    return names;
}

Properties Properties::excludeCopyProperties(
    std::vector<std::size_t> const& exclude_elem_ids,
    std::vector<std::size_t> const& exclude_node_ids) const
{
    Properties exclude_copy;
    for (auto name_vector_pair : properties_)
    {
        if (name_vector_pair.second->getMeshItemType() == MeshItemType::Cell)
        {
            exclude_copy.properties_.insert(std::make_pair(
                name_vector_pair.first,
                name_vector_pair.second->clone(exclude_elem_ids)));
        }
        else if (name_vector_pair.second->getMeshItemType() ==
                 MeshItemType::Node)
        {
            exclude_copy.properties_.insert(std::make_pair(
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
    for (auto name_vector_pair : properties_)
    {
        if (std::find(exclude_mesh_item_types.begin(),
                      exclude_mesh_item_types.end(),
                      name_vector_pair.second->getMeshItemType()) !=
            exclude_mesh_item_types.end())
        {
            continue;
        }

        std::vector<std::size_t> const exclude_positions{};
        new_properties.properties_.insert(
            std::make_pair(name_vector_pair.first,
                           name_vector_pair.second->clone(exclude_positions)));
    }
    return new_properties;
}

Properties::Properties(Properties const& properties)
    : properties_(properties.properties_)
{
    std::vector<std::size_t> exclude_positions;
    for (auto& name_vector_pair : properties_)
    {
        PropertyVectorBase* t(
            name_vector_pair.second->clone(exclude_positions));
        name_vector_pair.second = t;
    }
}

Properties& Properties::operator=(Properties const& properties)
{
    if (&properties == this)
    {
        return *this;
    }

    properties_ = properties.properties_;
    std::vector<std::size_t> exclude_positions;
    for (auto& name_vector_pair : properties_)
    {
        PropertyVectorBase* t(
            name_vector_pair.second->clone(exclude_positions));
        name_vector_pair.second = t;
    }

    return *this;
}

Properties::~Properties()
{
    for (auto name_vector_pair : properties_)
    {
        delete name_vector_pair.second;
    }
}

} // end namespace MeshLib

