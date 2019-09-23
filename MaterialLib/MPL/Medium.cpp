/**
 * \file
 * \author Norbert Grunwald
 * \date   07.09.2017
 *
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "Medium.h"

#include "BaseLib/Algorithm.h"
#include "Properties/Properties.h"

namespace MaterialPropertyLib
{
Medium::Medium(std::vector<std::unique_ptr<Phase>>&& phases,
               std::unique_ptr<PropertyArray>&& properties)
    : _phases(std::move(phases))
{
    if (properties)
    {
        overwriteExistingProperties(_properties, *properties, this);
    }
}

Phase const& Medium::phase(std::size_t const index) const
{
    return *_phases[index];
}

Phase const& Medium::phase(std::string const& name) const
{
    return *BaseLib::findElementOrError(
        _phases.begin(), _phases.end(),
        [&name](std::unique_ptr<MaterialPropertyLib::Phase> const& phase) {
            return phase->name() == name;
        },
        "Could not find phase name '" + name + "'.");
}

Property const& Medium::property(PropertyType const& p) const
{
    return *_properties[p];
}

std::size_t Medium::numberOfPhases() const
{
    return _phases.size();
}
}  // namespace MaterialPropertyLib
