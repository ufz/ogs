// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <memory>
#include <span>
#include <string>
#include <vector>

#include "Component.h"

namespace MaterialPropertyLib
{
class Property;

/// Enumeration of phase types.
enum class PhaseName
{
    Solid,
    AqueousLiquid,
    Gas,
    FrozenLiquid
};

/// Convert phase enum to its string representation.
[[nodiscard]] std::string_view toString(PhaseName phase_name);

/// Convert string to phase enum. Throws if invalid phase name.
[[nodiscard]] PhaseName fromString(std::string const& phase_name);

/// This class defines material phases.
///
/// The Phase class consists of a vector of components and an array of
/// properties.
class Phase final
{
public:
    /// The Phase constructor is called with the phase type enum.
    Phase(PhaseName phase_name,
          std::vector<std::unique_ptr<Component>>&& components,
          std::unique_ptr<PropertyArray>&& properties);

    /// A simple get-function for a component. The argument refers to the
    /// Index of the component in the components vector.
    Component const& component(std::size_t const& index) const;
    bool hasComponent(std::size_t const& index) const;

    /// A get-function for a component by component name.
    Component const& component(std::string const& name) const;

    /// A get-function for a property. The argument refers to the name of the
    /// property.
    Property const& property(PropertyType const& p) const;

    Property const& operator[](PropertyType const& p) const;

    bool hasProperty(PropertyType const& p) const;

    /// A get-function for retrieving the number of components in this phase.
    std::size_t numberOfComponents() const;

    /// Short description of the phase with its name.
    std::string description() const;

public:
    PhaseName const phaseName;

private:
    std::vector<std::unique_ptr<Component>> const components_;

    /// Here, all for all properties listed in the Properties enumerator are
    /// initialized by mole average functions of value zero. However,
    /// 'special-default' properties are allowed to be set.
    ///
    /// After this, other special properties can be set as exceptional defaults.
    PropertyArray properties_;
};

void checkRequiredProperties(
    Phase const& phase,
    std::span<PropertyType const> const required_properties);

}  // namespace MaterialPropertyLib
