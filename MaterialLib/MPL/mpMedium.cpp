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

namespace MaterialPropertyLib
{

Medium::Medium(BaseLib::ConfigTree const& config)
{
    // A Medium consists of phases and properties only.
    // Parse the phase configurations and push them into the
    // Medium::_phases attribute;
    auto const phases_config = config.getConfigSubtree("phases");
    createPhases(phases_config);

    // Parse the property configurations. These are medium properties only.
    // The properties of the phases are handled in the respective phases
    auto const properties_config = config.getConfigSubtree("properties");
    createProperties(properties_config);
}

void Medium::createPhases(BaseLib::ConfigTree const& config)
{
    std::vector<Phase*> phases;
    for (auto phase_config : config.getConfigSubtreeList("phase"))
    {
        // Unlike a medium, a phase may have a name. However, this is
        // silly at the moment since this name has no effect (except of some
        // benefits in terms of readability)
        auto const phase_name = phase_config.getConfigParameterOptional<std::string>("name");
//        Phase newPhase (phase_name);
        Phase* newPhase = new Phase(phase_name);
        // Furthermore, a phase (similar to a medium) consists of components and
        // properties.
        // Parsing the components:
        auto const components_config = phase_config.getConfigSubtree("components");
        newPhase->createComponents (components_config);
        // Parsing the phase properties:
        auto const properties_config = phase_config.getConfigSubtree("properties");
        newPhase->createProperties (properties_config);
        phases.push_back(newPhase);
    }
    _phases = phases;
}

void Medium::createProperties(BaseLib::ConfigTree const& config)
{
    for (auto property_config : config.getConfigSubtreeList("property"))
    {
        /// create a new Property based on configuration tree
        Property* property = newProperty (property_config);
        /// parse the name of the property
        auto const property_name =
                property_config.getConfigParameter<std::string>("name");
        /// insert the newly created property at the right place
        /// into the property array
        properties[convertStringToProperty(property_name)];
        _properties[convertStringToProperty(property_name)]=property;
    }
}

Phase* Medium::phase(std::size_t const index)
{
    return _phases[index];
}
std::size_t Medium::numberOfPhases(void)
{
    return _phases.size();
}

void Medium::summary()
{
	auto const nPhase = numberOfPhases();
	std::cout << "Number of phases: " << nPhase << "\n";
	for (size_t p=0; p < nPhase; ++p)
	{
		std::cout << "   Phase number " << p << ":\n";
		auto const nComponents = _phases[p]->numberOfComponents();
		std::cout << "   Number of Components: " << nComponents << "\n";
		for (size_t c=0; c < nComponents; ++c)
		{
			std::cout << "      Component number " << c << ":\n";
			std::string component_name =
					boost::get<std::string>
					(_phases[p]->component(c)->property(name)->value());
			std::cout << component_name << "\n";
		}
	}
}

} // MaterialPropertyLib




