/**
 * \file
 * \author Norbert Grunwald
 * \date   Sep 7, 2017
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "mpComponent.h"
#include "Components/components.h"
#include "Properties/properties.h"

namespace MaterialPropertyLib
{
Component::Component()
{
    // Default properties are set to zero.
    createDefaultProperties();

    // Some properties can be initialized by other default properties:
    _properties[name] = std::make_unique<Constant>("no_name");
}
Component::Component(std::string const& component_name)
{
    // Default properties are set to zero.
    createDefaultProperties();

    // Some properties can be initialized by other default properties:
    _properties[name] = std::make_unique<Constant>(component_name);
}

void Component::createProperties(BaseLib::ConfigTree const& config)
{
    for (auto property_config : config.getConfigSubtreeList("property"))
    {
        // Parsing the property name:
        auto const property_name =
            property_config.getConfigParameter<std::string>("name");
        // Create a new property based on the configuration subtree:
        auto property = newProperty(property_config, this);

        // Insert the new property at the right position into the components
        // private PropertyArray:
        _properties[convertStringToProperty(property_name)] =
            std::move(property);
    }
}

void Component::createDefaultProperties()
{
    for (std::size_t i = 0; i < number_of_property_enums; ++i)
    {
        _properties[i] =
            std::make_unique<Undefined>(static_cast<PropertyEnum>(i));
    }
}

Property& Component::property(PropertyEnum const& p) const
{
    return *_properties[p];
}

std::unique_ptr<Component> newComponent(std::string const& component_name,
                                        bool& isCustomComponent)
{
    // If a name is given, it must conform with one of the derived component
    // names in the following list:
    if (boost::iequals(component_name, "water"))
    {
        return std::make_unique<Water>();
    }

    // If the given name does not appear in the list (which is common), the
    // method creates a custom component, where all properties have to be
    // specified manually.
    isCustomComponent = true;
    return std::make_unique<Component>(component_name);
}
}  // namespace MaterialPropertyLib
