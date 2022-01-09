/**
 * \file
 * \author Norbert Grunwald
 * \date   Sep 7, 2017
 *
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "Component.h"

#include "Components/Components.h"
#include "Properties/Properties.h"

namespace MaterialPropertyLib
{
Component::Component() {}

Component::Component(std::string const& component_name,
                     std::unique_ptr<PropertyArray>&& properties)
    : name(component_name)
{
    if (properties)
    {
        overwriteExistingProperties(properties_, *properties, this);
    }
}

Property const& Component::property(PropertyType const& p) const
{
    Property const* const property = properties_[p].get();
    if (property == nullptr)
    {
        OGS_FATAL("Trying to access undefined property '{:s}' of {:s}",
                  property_enum_to_string[p], description());
    }

    return *properties_[p];
}

Property const& Component::operator[](PropertyType const& p) const
{
    return property(p);
}

bool Component::hasProperty(PropertyType const& p) const
{
    return properties_[p] != nullptr;
}

std::string Component::description() const
{
    return "component '" + name + "'";
}
}  // namespace MaterialPropertyLib
