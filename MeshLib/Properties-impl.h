/**
 * \file
 * \brief  Implemenatiom of the template part of the class Properties.
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

template <typename T>
PropertyVector<T>* Properties::createNewPropertyVector(
    std::string const& name,
    MeshItemType mesh_item_type,
    std::size_t n_components)
{
    std::map<std::string, PropertyVectorBase*>::const_iterator it(
        _properties.find(name)
    );
    if (it != _properties.end()) {
        ERR("A property of the name '{:s}' is already assigned to the mesh.",
            name);
        return nullptr;
    }
    auto entry_info(
        _properties.insert(
            std::make_pair(
                name, new PropertyVector<T>(name, mesh_item_type, n_components)
            )
        )
    );
    return static_cast<PropertyVector<T>*>((entry_info.first)->second);
}

template <typename T>
PropertyVector<T>* Properties::createNewPropertyVector(
    std::string const& name,
    std::size_t n_prop_groups,
    std::vector<std::size_t> const& item2group_mapping,
    MeshItemType mesh_item_type,
    std::size_t n_components)
{
    // check if there is already a PropertyVector with the same name
    std::map<std::string, PropertyVectorBase*>::const_iterator it(
        _properties.find(name)
    );
    if (it != _properties.end()) {
        ERR("A property of the name '{:s}' already assigned to the mesh.",
            name);
        return nullptr;
    }

    // check entries of item2group_mapping for consistence
    for (std::size_t k(0); k<item2group_mapping.size(); k++) {
        std::size_t const group_id (item2group_mapping[k]);
        if (group_id >= n_prop_groups) {
            ERR("The mapping to property {:d} for item {:d} is not in the "
                "correct range [0,{:d}).",
                group_id, k, n_prop_groups);
            return nullptr;
        }
    }

    auto entry_info(
        _properties.insert(
            std::pair<std::string, PropertyVectorBase*>(
                name,
                new PropertyVector<T>(n_prop_groups,
                    item2group_mapping, name, mesh_item_type, n_components)
            )
        )
    );
    return static_cast<PropertyVector<T>*>((entry_info.first)->second);
}

template <typename T>
bool Properties::existsPropertyVector(std::string const& name) const
{
    auto it(_properties.find(name));
    // Check that a PropertyVector with the appropriate name exists.
    if (it == _properties.end())
    {
        return false;
    }
    // Check that the PropertyVector has the correct data type.
    return dynamic_cast<PropertyVector<T> const*>(it->second) != nullptr;
}

template <typename T>
bool Properties::existsPropertyVector(std::string const& name,
                                      MeshItemType const mesh_item_type,
                                      int const number_of_components) const
{
    auto const it = _properties.find(name);
    if (it == _properties.end())
    {
        return false;
    }

    auto property = dynamic_cast<PropertyVector<T>*>(it->second);
    if (property == nullptr)
    {
        return false;
    }
    if (property->getMeshItemType() != mesh_item_type)
    {
        return false;
    }
    if (property->getNumberOfGlobalComponents() != number_of_components)
    {
        return false;
    }
    return true;
}

template <typename T>
PropertyVector<T> const* Properties::getPropertyVector(
    std::string const& name) const
{
    auto it(_properties.find(name));
    if (it == _properties.end())
    {
        OGS_FATAL("The PropertyVector '{:s}' is not available in the mesh.",
                  name);
    }
    if (!dynamic_cast<PropertyVector<T> const*>(it->second))
    {
        OGS_FATAL(
            "The PropertyVector '{:s}' has a different type than the requested "
            "PropertyVector.",
            name);
    }
    return dynamic_cast<PropertyVector<T> const*>(it->second);
}

template <typename T>
PropertyVector<T>* Properties::getPropertyVector(std::string const& name)
{
    auto it(_properties.find(name));
    if (it == _properties.end())
    {
        OGS_FATAL(
            "A PropertyVector with the specified name '{:s}' is not available.",
            name);
    }
    if (!dynamic_cast<PropertyVector<T>*>(it->second))
    {
        OGS_FATAL(
            "The PropertyVector '{:s}' has a different type than the requested "
            "PropertyVector.",
            name);
    }
    return dynamic_cast<PropertyVector<T>*>(it->second);
}

template <typename T>
PropertyVector<T> const* Properties::getPropertyVector(
    std::string const& name, MeshItemType const item_type,
    int const n_components) const
{
    auto const it = _properties.find(name);
    if (it == _properties.end())
    {
        OGS_FATAL(
            "A PropertyVector with name '{:s}' does not exist in the mesh.",
            name);
    }

    auto property = dynamic_cast<PropertyVector<T>*>(it->second);
    if (property == nullptr)
    {
        OGS_FATAL(
            "Could not cast the data type of the PropertyVector '{:s}' to "
            "requested data type.",
            name);
    }
    if (property->getMeshItemType() != item_type)
    {
        OGS_FATAL(
            "The PropertyVector '{:s}' has type '{:s}'. A '{:s}' field is "
            "requested.",
            name, toString(property->getMeshItemType()), toString(item_type));
    }
    if (property->getNumberOfGlobalComponents() != n_components)
    {
        OGS_FATAL(
            "PropertyVector '{:s}' has {:d} components, {:d} components are "
            "needed.",
            name, property->getNumberOfGlobalComponents(), n_components);
    }
    return property;
}

template <typename T>
PropertyVector<T>* Properties::getPropertyVector(std::string const& name,
                                                 MeshItemType const item_type,
                                                 int const n_components)
{
    auto const it = _properties.find(name);
    if (it == _properties.end())
    {
        OGS_FATAL(
            "A PropertyVector with name '{:s}' does not exist in the mesh.",
            name);
    }

    auto property = dynamic_cast<PropertyVector<T>*>(it->second);
    if (property == nullptr)
    {
        OGS_FATAL(
            "Could not cast the data type of the PropertyVector '{:s}' to "
            "requested data type.",
            name);
    }
    if (property->getMeshItemType() != item_type)
    {
        OGS_FATAL(
            "The PropertyVector '{:s}' has type '{:s}'. A '{:s}' field is "
            "requested.",
            name, toString(property->getMeshItemType()), toString(item_type));
    }
    if (property->getNumberOfGlobalComponents() != n_components)
    {
        OGS_FATAL(
            "PropertyVector '{:s}' has {:d} components, {:d} components are "
            "needed.",
            name, property->getNumberOfGlobalComponents(), n_components);
    }
    return property;
}
