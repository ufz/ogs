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

namespace MaterialPropertyLib
{

Phase::Phase(){};
Phase::Phase(boost::optional<std::string> const& name)
{
    if (name)
        _name = name.get();
    else
        _name = "unknown";
};

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
    {
    }
}


} // MaterialPropertyLib

