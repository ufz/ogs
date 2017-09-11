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

#include "mpComponent.h"
#include "Components/components.h"


namespace MaterialPropertyLib
{

Component::Component(){};

void Component::createProperties(BaseLib::ConfigTree const& config)
{
    for (auto property_config : config.getConfigSubtreeList("property"))
    {
        /// Properties always consist of a name and a type. The name in this
        /// case refers to the sort of property (such as density, molar_mass,
        /// etc.),while the type identifies the specific function which is
        /// used to compute the property. "type" therefore always refers to
        /// a specific derived class of the base class Property (such as
        /// IdealGasLaw, Constant, PengRobinson, etc.).

        /// Create a new property based on the configuration subtree:
        Property* property = newProperty (property_config);

        /// Parsing the property name:
        auto const property_name =
                property_config.getConfigParameter<std::string>("name");

        /// Insert the new property at the right position into the components
        /// private PropertyArray:
        _properties[convertStringToProperty(property_name)] = property;
    }
}

Property* Component::property(PropertyEnum const &p)
{
    return _properties[p];
}

Component* newComponent (boost::optional<std::string> const &name)
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


