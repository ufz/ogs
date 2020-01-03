/**
 * \file
 * \author Norbert Grunwald
 * \date   Sep 7, 2017
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "Component.h"
#include "Components/Components.h"
#include "Properties/Properties.h"

namespace MaterialPropertyLib
{
Component::Component()
{
    // Some properties can be initialized by other default properties:
    _properties[PropertyType::name] = std::make_unique<Constant>("no_name");
}

Component::Component(std::string const& component_name,
                     std::unique_ptr<PropertyArray>&& properties)
{
    // Some properties can be initialized by other default properties:
    _properties[PropertyType::name] =
        std::make_unique<Constant>(component_name);

    if (properties)
    {
        overwriteExistingProperties(_properties, *properties, this);
    }
}

Property const& Component::property(PropertyType const& p) const
{
    return *_properties[p];
}

std::string Component::name() const
{
    return std::get<std::string>(_properties[PropertyType::name]->value());
}
}  // namespace MaterialPropertyLib
