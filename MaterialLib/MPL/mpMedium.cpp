/**
 * \author Norbert Grunwald
 * \date   07.09.2017
 * \brief
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "mpMedium.h"
#include "mpPhase.h"
#include "mpComponent.h"
#include "Properties/pAverageVolumeFraction.h"
#include "Properties/pConstant.h"
#include <iostream>
#include <string>
#include <boost/variant.hpp>

namespace MaterialPropertyLib
{
/**
 * This constructor parses the "phases" and "properties" subtrees
 * of the config tree and calls create methods for the phase vector
 * and the properties array. Medium properties are optional. If not
 * defined, default properties are assigned.
 */
Medium::Medium(BaseLib::ConfigTree const& config)
{
    // Default properties are initialized in the first place, such
    // that user-defined properties overwrite those defaults.
    createDefaultProperties();

    // A name of a medium is entirely optional and has no effects
    // other than a small gain of clarity in case several media
    // should be defined.
    auto const medium_name =
            config.getConfigParameterOptional<std::string>("name");
    if (medium_name)
            _properties[name] = new Constant(medium_name.get());
        else
            _properties[name] = new Constant("no_name");

	// Parsing the phases
    auto const phases_config = config.getConfigSubtree("phases");
    createPhases(phases_config);
    // Parsing medium properties, overwriting the defaults.
    auto const properties_config =
            config.getConfigSubtreeOptional("properties");
    if (properties_config)
    	createProperties(properties_config.get());
}
/**
 * This method creates the phases of the medium. Unlike a medium, a
 * phase may have a name. However, this is silly at the moment since
 * this name still has no effect (except of some benefits in regard of
 * readability).
 * Phase components are required (a phase consists of at least one
 * component).
 * Phase properties are optional. If not given, default properties
 * are assigned. These default properties average the component
 * properties, weighted by mole fraction.
 */
void Medium::createPhases(BaseLib::ConfigTree const& config)
{
    std::vector<Phase*> phases;
    for (auto phase_config : config.getConfigSubtreeList("phase"))
    {
        // Phase name is optional
        auto const phase_name = phase_config.getConfigParameterOptional<std::string>("name");
        Phase* newPhase = new Phase(phase_name);
        // Parsing the components
        auto const components_config = phase_config.getConfigSubtree("components");
        newPhase->createComponents (components_config);

        // Properties of phases are optional
        if (auto const properties_config = phase_config.getConfigSubtreeOptional("properties"))
                	newPhase->createProperties (properties_config.get());
        // No else branch here, default properties are used. Those defaults
        // were assigned by the phase constructor.
        phases.push_back(newPhase);
    }
    _phases = phases;
}
/**
 * This method creates the properties of the Medium as defined in the
 * prj-file. Only specified properties overwrite the default properties.
 */
void Medium::createProperties(BaseLib::ConfigTree const& config)
{
    for (auto property_config : config.getConfigSubtreeList("property"))
    {
        // create a new Property based on configuration tree
        Property* property = newProperty (property_config, this);
        /// parse the name of the property
        auto const property_name =
                property_config.getConfigParameter<std::string>("name");
        // insert the newly created property at the right place
        // into the property array
        _properties[convertStringToProperty(property_name)]=property;
    }
}

/**
 * This method defines default properties of the medium. Most properties
 * are fine with the volume fraction average, but special-defaults are
 * allowed as well...
 */
void Medium::createDefaultProperties(void)
{
	for (size_t i=0; i < number_of_property_enums; ++i)
		this->_properties[i] = new AverageVolumeFraction(this);
}

Phase* Medium::phase(std::size_t const index)
{
    return _phases[index];
}

Property* Medium::property(PropertyEnum const &p)
{
    return _properties[p];
}

std::size_t Medium::numberOfPhases(void)
{
    return _phases.size();
}

void Medium::resetPropertyValues()
{
    for (size_t i=0; i<number_of_property_enums; ++i)
    {
        _properties[i]->set(false);
    }
}
void Medium::summary()
{
    VariableArray stdRefCond;

    stdRefCond[T] = 293.15;
    stdRefCond[p_GR] = 101325.;
    stdRefCond[p_LR] = 101325.;

    std::cout << "\n";
    std::cout << "Material summary: Showing only non-zero properties!\n";
    std::cout << "==================================================\n";

    std::string mediumName = getString(this->property(name));
     std::cout << "Medium: \'" << mediumName << "\'\n";
    std::cout << "----------";
    for (short i=0; i<mediumName.length(); i++)
        std::cout << "-";
    std::cout << "\n";

	auto const nPhase = numberOfPhases();
	std::cout << "Number of phases: " << nPhase << "\n";

	for (size_t p=0; p < nPhase; ++p)
	{
	    Phase* thisPhase = _phases[p];
	    std::string phaseName = getString(thisPhase->property(name));
		std::cout << "   Phase number " << p << ": \'" << phaseName << "\'\n";
	    std::cout << "   -----------------";
	    for (short i=0; i<phaseName.length(); i++)
	        std::cout << "-";
	    std::cout << "\n";
		auto const nComponents = _phases[p]->numberOfComponents();
		std::cout << "   Number of Components: " << nComponents << "\n";
		for (size_t c=0; c < nComponents; ++c)
		{
		    Component* thisComponent = thisPhase->component(c);
		    std::string componentName =
		            getString(thisComponent->property(name));

		    std::cout << "      Component number " << c << ": \'" << componentName << "\'\n";
	        std::cout << "      ---------------------";
	        for (short i=0; i<componentName.length(); i++)
	            std::cout << "-";
	        std::cout << "\n";
		    for (size_t prop = 0; prop < number_of_property_enums; ++prop)
		    {
		        Property* thisProperty = thisComponent->property(static_cast<PropertyEnum>(prop));
		        auto property_name = convertEnumToString[prop];
		        if (property_name == "name")
		            continue;
		        const double property_value = getScalar(thisProperty);
		        if (property_value != 0)
		        std::cout << "         " << property_name << ": "
		                << property_value << "\n";

		    }
		}
        for (size_t prop = 0; prop < number_of_property_enums; ++prop)
        {
            Property* thisProperty = thisPhase->property(static_cast<PropertyEnum>(prop));
            auto property_name = convertEnumToString[prop];
            if (property_name == "name")
                continue;
            const double property_value = getScalar(thisProperty);
            if (property_value != 0)
            std::cout << "      " << property_name << ": "
                    << property_value << "\n";
        }
	}
    for (size_t prop = 0; prop < number_of_property_enums; ++prop)
    {
        Property* thisProperty = property(static_cast<PropertyEnum>(prop));
        auto property_name = convertEnumToString[prop];
        if (property_name == "name")
            continue;
        const double property_value = getScalar(thisProperty);
        if (property_value != 0)
        std::cout << "   " << property_name << ": "
                << property_value << "\n";
    }

}

} // MaterialPropertyLib




