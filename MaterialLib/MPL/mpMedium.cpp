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
            config.getConfigParameter<std::string>("name","no_name");

    _properties[name] = std::make_unique<Constant>(medium_name);

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
    for (auto phase_config : config.getConfigSubtreeList("phase"))
    {
        // Phase name is optional
        auto const phase_name = phase_config.getConfigParameterOptional<std::string>("name");
        auto newPhase = std::make_unique<Phase>(phase_name);
        // Parsing the components
        auto const components_config = phase_config.getConfigSubtree("components");
        newPhase->createComponents (components_config);

        // Properties of phases are optional
        if (auto const properties_config = phase_config.getConfigSubtreeOptional("properties"))
                	newPhase->createProperties (properties_config.get());
        // No else branch here, default properties are used. Those defaults
        // were assigned by the phase constructor.
        _phases.push_back(std::move(newPhase));
    }
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
        auto property = newProperty (property_config, this);
        /// parse the name of the property
        auto const property_name =
                property_config.getConfigParameter<std::string>("name");
        // insert the newly created property at the right place
        // into the property array
        _properties[convertStringToProperty(property_name)]= std::move(property);
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
		this->_properties[i] = std::make_unique<AverageVolumeFraction>(this);
}

Phase* Medium::phase(std::size_t const index)
{
    return _phases[index].get();
}

Property* Medium::property(PropertyEnum const &p)
{
    return _properties[p].get();
}

std::size_t Medium::numberOfPhases(void)
{
    return _phases.size();
}

void Medium::resetPropertyUpdateStatus()
{
    // Phase properties
    for (size_t p = 0; p < _phases.size(); ++p)
        _phases[p]->resetPropertyUpdateStatus();
    // Medium properties
    for (size_t i=0; i<number_of_property_enums; ++i)
        _properties[i]->isUpdated(false);
}

} // MaterialPropertyLib




