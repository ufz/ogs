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

{
        {
        }
}

} // MaterialPropertyLib


