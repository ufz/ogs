/**
 * \file
 * \author Norbert Grunwald
 * \date   Sep 7, 2017
 *
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "Phase.h"
#include "BaseLib/Algorithm.h"
#include "Properties/Properties.h"

#include "Component.h"

namespace MaterialPropertyLib
{
Phase::Phase(std::string&& phase_name,
             std::vector<std::unique_ptr<Component>>&& components,
             std::unique_ptr<PropertyArray>&& properties)
    : _components(std::move(components))
{
    _properties[PropertyType::name] = std::make_unique<Constant>(phase_name);

    if (properties)
    {
        overwriteExistingProperties(_properties, *properties, this);
    }
}

Component const& Phase::component(const std::size_t& index) const
{
    return *_components[index];
}

Component const& Phase::component(std::string const& name) const
{
    return *BaseLib::findElementOrError(
        _components.begin(), _components.end(),
        [&name](std::unique_ptr<Component> const& component) {
            return component->name() == name;
        },
        "Could not find component name '" + name + "'.");
}

Property const& Phase::property(PropertyType const& p) const
{
    return *_properties[p];
}

std::size_t Phase::numberOfComponents() const
{
    return _components.size();
}

std::string Phase::name() const
{
    return std::get<std::string>(_properties[PropertyType::name]->value());
}
}  // namespace MaterialPropertyLib
