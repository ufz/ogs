/**
 * \file
 * \author Norbert Grunwald
 * \date   Sep 7, 2017
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "Phase.h"

#include "BaseLib/Algorithm.h"
#include "Component.h"
#include "Properties/Properties.h"

namespace MaterialPropertyLib
{
Phase::Phase(std::string&& phase_name,
             std::vector<std::unique_ptr<Component>>&& components,
             std::unique_ptr<PropertyArray>&& properties)
    : name(std::move(phase_name)), components_(std::move(components))
{
    if (properties)
    {
        overwriteExistingProperties(properties_, *properties, this);
    }
}

Component const& Phase::component(const std::size_t& index) const
{
    return *components_[index];
}

bool Phase::hasComponent(std::size_t const& index) const
{
    return components_[index] != nullptr;
}

Component const& Phase::component(std::string const& name) const
{
    return *BaseLib::findElementOrError(
        components_.begin(), components_.end(),
        [&name](std::unique_ptr<Component> const& component) {
            return component->name == name;
        },
        "Could not find component name '" + name + "'.");
}

Property const& Phase::property(PropertyType const& p) const
{
    Property const* const property = properties_[p].get();
    if (property == nullptr)
    {
        OGS_FATAL("Trying to access undefined property '{:s}' of {:s}",
                  property_enum_to_string[p], description());
    }
    return *properties_[p];
}

Property const& Phase::operator[](PropertyType const& p) const
{
    return property(p);
}

bool Phase::hasProperty(PropertyType const& p) const
{
    return properties_[p] != nullptr;
}

std::size_t Phase::numberOfComponents() const
{
    return components_.size();
}

std::string Phase::description() const
{
    return "phase '" + name + "'";
}
}  // namespace MaterialPropertyLib
