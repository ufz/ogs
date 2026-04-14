// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "Phase.h"

#include "BaseLib/Algorithm.h"
#include "BaseLib/Error.h"
#include "Component.h"
#include "Properties/Properties.h"

namespace MaterialPropertyLib
{
std::string_view toString(PhaseName phase_name)
{
    switch (phase_name)
    {
        case PhaseName::Solid:
            return "Solid";
        case PhaseName::AqueousLiquid:
            return "AqueousLiquid";
        case PhaseName::Gas:
            return "Gas";
        case PhaseName::FrozenLiquid:
            return "FrozenLiquid";
    }
    OGS_FATAL("Unknown phase name");
}

PhaseName fromString(std::string const& phase_name)
{
    if (phase_name == "Solid")
    {
        return PhaseName::Solid;
    }
    if (phase_name == "AqueousLiquid")
    {
        return PhaseName::AqueousLiquid;
    }
    if (phase_name == "Gas")
    {
        return PhaseName::Gas;
    }
    if (phase_name == "FrozenLiquid")
    {
        return PhaseName::FrozenLiquid;
    }
    OGS_FATAL("Unknown phase name '{}'", phase_name);
}

Phase::Phase(PhaseName phase_name,
             std::vector<std::unique_ptr<Component>>&& components,
             std::unique_ptr<PropertyArray>&& properties)
    : phaseName(phase_name), components_(std::move(components))
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
        components_,
        [&name](
            std::unique_ptr<MaterialPropertyLib::Component> const& component)
        { return component->name == name; },
        [&]() { OGS_FATAL("Could not find component named '{:s}'.", name); });
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
    return "phase '" + std::string(toString(phaseName)) + "'";
}

void checkRequiredProperties(
    Phase const& phase, std::span<PropertyType const> const required_properties)
{
    for (auto const& p : required_properties)
    {
        if (!phase.hasProperty(p))
        {
            OGS_FATAL("The property '{:s}' is missing in the {:s} phase.",
                      property_enum_to_string[p], toString(phase.phaseName));
        }
    }
}

}  // namespace MaterialPropertyLib
