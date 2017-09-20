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
#include "Properties/properties.h"


namespace MaterialPropertyLib
{
/// Default constructor of Component. This constructor is used
/// when the component is not specified via the 'name'-tag.
Component::Component()
{
	DBUG ("      Unspecified component constructed.");
	// Default properties are set to zero.
	createDefaultProperties();

	// Some properties can be initialized by other default
	// properties:
	_properties[name] = new Constant ("no_name");
}
/**
 * \brief The method reads the 'properties' tag in the prj-file and creates
 * component properties accordingly.
 * \details First, a new property iy created based on the specified property
 * type. Then, the property name is evaluated and the property is copied
 * into the properties array.
 */
void Component::createProperties(BaseLib::ConfigTree const& config)
{
    for (auto property_config : config.getConfigSubtreeList("property"))
    {
        /// Create a new property based on the configuration subtree:
        Property* property = newProperty (property_config, this);

        /// Parsing the property name:
        auto const property_name =
                property_config.getConfigParameter<std::string>("name");

        /// Insert the new property at the right position into the components
        /// private PropertyArray:
        _properties[convertStringToProperty(property_name)] = property;
    }
}

/**
 * This method initializes the component property array with
 * constant functions of value zero.
 */
void Component::createDefaultProperties(void)
{
for (size_t i=0; i < number_of_property_enums; ++i)
	this->_properties[i] = new Constant(0.);
}

Property* Component::property(PropertyEnum const &p)
{
    return _properties[p];
}
/**
 * \brief This function creates a new component based on the (optional)
 * component name that is given in the prj-file.
 * \details The method evaluates the string in the 'name'-object and
 * calls the constructors of the derived component classes (if found)
 * or that of the base class (if no name is specified).
 */
Component* newComponent (boost::optional<std::string> const &name)
{
    // Check whether a name is given or not
	if (name)
	{
		std::string component_name = static_cast<std::string>(name.get());
		// If a name is given, it must conform with one of the
		// derived component names in the following list:
		if (boost::iequals(component_name, "water")) // "string" == "StRiNg"
		{
			return new Water;
		}
		if (boost::iequals(component_name, "salt"))
		{
			return new Salt;
		}
		if (boost::iequals(component_name, "carbondioxide"))
		{
			return new CarbonDioxide;
		}
		else
		{
	        // If the given name does not appear in the list, this
	        // throws a fatal error.
			OGS_FATAL("The specified component name \"%s\" was not recognized", component_name.c_str());
			return nullptr; // to avoid the 'no return' warning.
		}
	}
    // If there is no component name given (which is common), the
    // method creates a custom component, where all properties have
    // to be specified manually.
    else return new Component;

}

} // MaterialPropertyLib


