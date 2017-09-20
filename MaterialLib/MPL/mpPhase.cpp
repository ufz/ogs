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
#include "mpComponent.h"
#include "Properties/properties.h"

namespace MaterialPropertyLib
{
Phase::Phase(boost::optional<std::string> const& phase_name)
{
   createDefaultProperties();

    if (phase_name)
    	_properties[name] = new Constant(phase_name.get());
    else
    	_properties[name] = new Constant("unknown");
};
Phase::Phase()
{
   createDefaultProperties();
  _properties[name] = new Constant("no_name");
};
/**
 *  Just like a phase, a component can have a name. But, in this case,
 *  the name has an important task. If a name is given, a specific
 *  component class referring to that name with corresponding physical
 *  material properties is created.
 *  Assigning a name is optional; If no name is given, a custom
 *  component without predefined properties is created.
 */
void Phase::createComponents(BaseLib::ConfigTree const& config)
{
	std::vector<Component*> components;
	for (auto component_config : config.getConfigSubtreeList("component"))
	{

		// Parsing the optional component name
		auto const component_name =
				component_config.getConfigParameterOptional<std::string>("name");
		Component* component = newComponent (component_name);

		// Parsing component properties. If a component name is given,
		// properties are optional (since they may be predefined in the
		// class implementation).
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
			if (!component_name)
			{
				OGS_FATAL ("No Properties defined for unspecified component");
			}
		}

		components.push_back(component);
	}
	_components = components;
}
/**
 * This method creates the properties of the Phase as defined in the
 * prj-file. Only specified properties overwrite the default properties.
 */
void Phase::createProperties(BaseLib::ConfigTree const& config)
{
    for (auto property_config : config.getConfigSubtreeList("property"))
    {
        // create a new Property based on configuration tree
        Property* property = newProperty (property_config, this);
        // parse the name of the property
        auto const property_name =
                property_config.getConfigParameter<std::string>("name");
        // insert the newly created property at the right place
        // into the property array
        _properties[convertStringToProperty(property_name)]= property;
    }
}

/**
 * The default phase properties method.
 * Here, all for all properties listed in the Properties enumerator
 * are initialized by mole average functions of value zero. However,
 * 'special-default' properties are allowed to be set.
 */
void Phase::createDefaultProperties(void)
{
	for (size_t i=0; i < number_of_property_enums; ++i)
		this->_properties[i] = new AverageMoleFraction(this);
	// After this, other special properties can be set as default
}

Component* Phase::component(const std::size_t &index)
{
    return _components[index];
}

Property* Phase::property(PropertyEnum const &p)
{
    return _properties[p];
}

std::size_t Phase::numberOfComponents(void)
{
    return _components.size();
}

} // MaterialPropertyLib

