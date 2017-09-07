/**
 * \author Norbert Grunwald
 * \date   Sep 7, 2017
 * \brief
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "mpPhase.h"
#include "Components/components.h"

namespace MaterialPropertyLib
{

Phase::Phase(){};
Phase::Phase(std::string name)
: _name (name){};

void Phase::createComponents(BaseLib::ConfigTree const& config)
{
    std::vector<Component*> components;
    for (auto component_config : config.getConfigSubtreeList("component"))
    {
        /// Just like a phase, a component can have a name. But, in this case,
        /// the name has an important task. If a name is given, a specific
        /// component class referring to that name with corresponding physical
        /// material properties is created.
        /// Assigning a name is optional; If no name is given, a custom
        /// component without predefined properties is created.
        auto const component_name =
                component_config.getConfigParameterOptional<std::string>("name");
        Component* component = newComponent (component_name);

        // Parsing component properties:
        auto const properties_config = component_config.getConfigSubtree("properties");
        component->createProperties(properties_config);

        components.push_back(component);
    }
    _components = components;
}

void Phase::createProperties(BaseLib::ConfigTree const& config)
{

}

Component* Phase::newComponent (boost::optional<std::string> const &name)
{
    /// Check whether a name is given or not
    if (name)
    {
        std::string component_name = static_cast<std::string>(name.get());
        /// If a name is given, it must conform with one of the
        /// derived component names in the following list:
        if (boost::iequals(component_name, "water")) // "string" == "StRiNg"
        {
            return new Water;
        }

        /// If the given name does not appear in the list, this
        /// throws a fatal error.
        else
            {
            OGS_FATAL("The specified component name \"%s\" was not recognized", component_name.c_str());
            return nullptr; // to avoid the 'no return' warning.
            }
    }
    /// If there is no component name given (which is common), the
    /// method creates a custom component, where all properties have
    /// to be specified manually.
    else return new Component;

}


} // MaterialPropertyLib

