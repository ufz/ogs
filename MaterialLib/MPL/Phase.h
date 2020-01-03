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
#pragma once

#include <memory>
#include <string>
#include <vector>

#include "Component.h"

namespace MaterialPropertyLib
{
class Property;
}
namespace MaterialPropertyLib
{
/// This class defines material phases.
///
/// The Phase class consists of a vector of components and an array of
/// properties.
class Phase final
{
public:
    /// The Phase constructor is called with the optional phase name.
    Phase(std::string&& phase_name,
          std::vector<std::unique_ptr<Component>>&& components,
          std::unique_ptr<PropertyArray>&& properties);

    /// A simple get-function for a component. The argument refers to the
    /// Index of the component in the components vector.
    Component const& component(std::size_t const& index) const;

    /// A get-function for a component by component name.
    Component const& component(std::string const& name) const;

    /// A get-function for a property. The argument refers to the name of the
    /// property.
    Property const& property(PropertyType const& p) const;
    bool hasProperty(PropertyType const& p) const;

    /// A get-function for retrieving the number of components in this phase.
    std::size_t numberOfComponents() const;

    std::string name() const;

private:
    std::vector<std::unique_ptr<Component>> const _components;

    /// Here, all for all properties listed in the Properties enumerator are
    /// initialized by mole average functions of value zero. However,
    /// 'special-default' properties are allowed to be set.
    ///
    /// After this, other special properties can be set as exceptional defaults.
    PropertyArray _properties;
};

}  // namespace MaterialPropertyLib
