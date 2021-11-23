/**
 * \file
 * \author Norbert Grunwald
 * \date   07.09.2017
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "Medium.h"

#include "BaseLib/Algorithm.h"
#include "Properties/Properties.h"

namespace MaterialPropertyLib
{
Medium::Medium(int const material_id,
               std::vector<std::unique_ptr<Phase>>&& phases,
               std::unique_ptr<PropertyArray>&& properties)
    : phases_(std::move(phases)), material_id_(material_id)
{
    if (properties)
    {
        overwriteExistingProperties(properties_, *properties, this);
    }
}

Phase const& Medium::phase(std::size_t const index) const
{
    return *phases_[index];
}

Phase const& Medium::phase(std::string const& phase_name) const
{
    return *BaseLib::findElementOrError(
        phases_.begin(), phases_.end(),
        [&phase_name](std::unique_ptr<MaterialPropertyLib::Phase> const& phase)
        { return phase->name == phase_name; },
        "Could not find phase name '" + phase_name + "'.");
}

bool Medium::hasPhase(std::string const& phase_name) const
{
    return std::any_of(begin(phases_), end(phases_),
                       [&phase_name](auto const& phase)
                       { return phase->name == phase_name; });
}

Property const& Medium::property(PropertyType const& p) const
{
    Property const* const property = properties_[p].get();
    if (property == nullptr)
    {
        OGS_FATAL("Trying to access undefined property '{:s}' of {:s}",
                  property_enum_to_string[p], description());
    }
    return *properties_[p];
}

Property const& Medium::operator[](PropertyType const& p) const
{
    return property(p);
}

bool Medium::hasProperty(PropertyType const& p) const
{
    return properties_[p] != nullptr;
}

std::size_t Medium::numberOfPhases() const
{
    return phases_.size();
}

std::string Medium::description() const
{
    return "medium " + std::to_string(material_id_);
}

Phase const& fluidPhase(Medium const& medium)
{
    if (medium.hasPhase("Gas"))
    {
        return medium.phase("Gas");
    }
    if (medium.hasPhase("AqueousLiquid"))
    {
        return medium.phase("AqueousLiquid");
    }
    OGS_FATAL(
        "Neither Gas nor AqueousLiquid phase is available for the medium, but "
        "a fluid phase was requested.");
}
}  // namespace MaterialPropertyLib
