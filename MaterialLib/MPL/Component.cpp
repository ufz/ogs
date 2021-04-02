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
    return *properties_[p];
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
