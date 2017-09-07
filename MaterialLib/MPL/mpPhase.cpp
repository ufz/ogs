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

#include "mpPhase.h"
#include "BaseLib/Algorithm.h"
#include "Properties/properties.h"

#include "mpComponent.h"

namespace MaterialPropertyLib
{
Phase::Phase(std::string const& phase_name)
{
    createDefaultProperties();

    std::array<std::string, 4> allowed_phase_names = {
        {"solid", "aqueous liquid", "non-aqueous liquid", "gas"}};

    if (std::find(allowed_phase_names.begin(),
                  allowed_phase_names.end(),
                  phase_name) == allowed_phase_names.end())
    {
        ERR("Phase name should be one of:");
        for (auto const name : allowed_phase_names)
        {
            ERR(name.c_str());
        }
        OGS_FATAL("Wrong phase name '%s' given.", phase_name.c_str());
    }
    _properties[name] = std::make_unique<Constant>(Constant(phase_name));
}

void Phase::createComponents(BaseLib::ConfigTree const& config)
{
    // collect names of components
    std::set<std::string> component_names;

    for (auto component_config : config.getConfigSubtreeList("component"))
    {
        // Parsing the component name
        auto const component_name =
            component_config.getConfigParameter<std::string>("name");

        // Check whether a name is given or not
        if (component_name.empty())
            OGS_FATAL(
                "Component name is a mandatory field and cannot be empty.");

        bool isCustomComponent = false;
        auto component = newComponent(component_name, isCustomComponent);

        // Parsing component properties. If a component name is given and this
        // component is predefined in the class implementation, properties
        // become optional. The default values of properties will be overwritten
        // if specified.
        if (auto const properties_config =
                component_config.getConfigSubtreeOptional("properties"))
        {
            component->createProperties(properties_config.get());
        }
        else
        {
            // No component properties are provided. If the component is
            // not specified, this results in a fatal error, since an
            // unspecified component has no properties.
            if (isCustomComponent)
            {
                OGS_FATAL("No Properties defined for unspecified component");
            }
        }

        component_names.insert(component_name);
        _components.push_back(std::move(component));
    }

    if (component_names.size() != _components.size())
    {
        OGS_FATAL(
            "Found duplicates with the same component name tag inside a "
            "phase.");
    }
}

void Phase::createProperties(BaseLib::ConfigTree const& config)
{
    for (auto property_config : config.getConfigSubtreeList("property"))
    {
        // create a new Property based on configuration tree
        auto property = newProperty(property_config, this);
        // parse the name of the property
        auto const property_name =
            property_config.getConfigParameter<std::string>("name");
        // insert the newly created property at the right place into the
        // property array
        _properties[convertStringToProperty(property_name)] =
            std::move(property);
    }
}

void Phase::createDefaultProperties()
{
    for (std::size_t i = 0; i < number_of_property_enums; ++i)
    {
        _properties[i] = std::make_unique<Undefined>((PropertyEnum)i);
    }

    // After this, other special properties can
    // be set as exceptional defaults
}

Component& Phase::component(const std::size_t& index) const
{
    return *_components[index];
}

Component& Phase::component(std::string const& name) const
{
    return *BaseLib::findElementOrError(
        _components.begin(), _components.end(),
        [&name](std::unique_ptr<Component> const& p) {
            return getString(p->property(PropertyEnum::name)) == name;
        },
        "Could not find component name '" + name + "'.");
}

Property& Phase::property(PropertyEnum const& p) const
{
    return *_properties[p];
}

std::size_t Phase::numberOfComponents() const
{
    return _components.size();
}
}  // namespace MaterialPropertyLib
