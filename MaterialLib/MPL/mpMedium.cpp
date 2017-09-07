/**
 * \file
 * \author Norbert Grunwald
 * \date   07.09.2017
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "mpMedium.h"

#include <boost/variant.hpp>
#include <string>

#include "BaseLib/Algorithm.h"

#include "Properties/pUndefined.h"
#include "Properties/pConstant.h"
#include "mpComponent.h"
#include "mpPhase.h"

namespace MaterialPropertyLib
{
Medium::Medium(BaseLib::ConfigTree const& config)
{
    // Default properties are initialized in the first place, such
    // that user-defined properties overwrite those defaults.
    createDefaultProperties();

    // A name of a medium is entirely optional and has no effects
    // other than a small gain of clarity in case several media
    // should be defined.
    auto const medium_name =
        config.getConfigParameter<std::string>("name", "no_name");

    // 'name' is a constant property of the medium
    _properties[name] = std::make_unique<Constant>(medium_name);

    // Parsing the phases
    // Properties of phases may be not required in all the cases.
    auto const phases_config = config.getConfigSubtreeOptional("phases");
    if (phases_config)
        createPhases(*phases_config);

    // Parsing medium properties, overwriting the defaults.
    auto const properties_config =
        config.getConfigSubtreeOptional("properties");
    if (properties_config)
    {
        createProperties(*properties_config);
    }

    if (!phases_config && !properties_config)
        OGS_FATAL("Neither tag <phases> nor tag <properties> has been found.");
}
void Medium::createPhases(BaseLib::ConfigTree const& config)
{
    std::set<std::string> phase_names;

    for (auto phase_config : config.getConfigSubtreeList("phase"))
    {
        auto const phase_name =
            phase_config.getConfigParameter<std::string>("name");

        if (phase_name.empty())
            OGS_FATAL("Phase name is a mandatory field and cannot be empty.");

        auto newPhase = std::make_unique<Phase>(phase_name);
        // Parsing the components
        auto const components_config =
            phase_config.getConfigSubtreeOptional("components");
        if (components_config)
            newPhase->createComponents(components_config.get());

        // Properties of phases are optional
        auto const properties_config =
            phase_config.getConfigSubtreeOptional("properties");
        if (properties_config)
            newPhase->createProperties(properties_config.get());

        if (!components_config && !properties_config)
            OGS_FATAL(
                "Neither tag <components> nor tag <properties> has been "
                "found.");

        phase_names.insert(phase_name);

        _phases.push_back(std::move(newPhase));
    }

    if (phase_names.size() != _phases.size())
        OGS_FATAL(
            "Found duplicates with the same phase name tag inside a medium.");
}
void Medium::createProperties(BaseLib::ConfigTree const& config)
{
    for (auto property_config : config.getConfigSubtreeList("property"))
    {
        // create a new Property based on configuration tree
        auto property = newProperty(property_config, this);
        /// parse the name of the property
        auto const property_name =
            property_config.getConfigParameter<std::string>("name");
        // insert the newly created property at the right place
        // into the property array
        _properties[convertStringToProperty(property_name)] =
            std::move(property);
    }
}

void Medium::createDefaultProperties()
{
    for (std::size_t i = 0; i < number_of_property_enums; ++i)
    {
        this->_properties[i] = std::make_unique<Undefined>((PropertyEnum)i);
    }
}

Phase& Medium::phase(std::size_t const index) const
{
    return *_phases[index];
}

Phase& Medium::phase(std::string const& name) const
{
    return *BaseLib::findElementOrError(
        _phases.begin(), _phases.end(),
        [&name](std::unique_ptr<MaterialPropertyLib::Phase> const& p) {
            return getString(p->property(
                       MaterialPropertyLib::PropertyEnum::name)) == name;
        },
        "Could not find phase name '" + name + "'.");
}

Property& Medium::property(PropertyEnum const& p) const
{
    return *_properties[p];
}

std::size_t Medium::numberOfPhases() const
{
    return _phases.size();
}

}  // namespace MaterialPropertyLib
