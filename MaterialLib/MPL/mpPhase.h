/**
 * \author Norbert Grunwald
 * \date   Sep 7, 2017
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */
#pragma once

#include <string>
#include <vector>

#include "BaseLib/ConfigTree.h"

#include "mpComponent.h"
#include "mpProperty.h"

namespace MaterialPropertyLib
{
/// This class defines material phases.
///
/// The Phase class consists of a vector of components and an array of
/// properties.
class Phase final
{
private:
    std::vector<std::unique_ptr<Component>> _components;
    PropertyArray _properties;

public:
    /// The Phase constructor is called with the optional phase name.
    explicit Phase(std::string const& phase_name);

    /// The method creating phase components based on config subtree.
    ///
    /// Just like a phase, a component can have a name. But, in this case, the
    /// name has an important task. If a name is given, a specific component
    /// class referring to that name with corresponding physical material
    /// properties is created.
    /// Assigning a name is optional; If no name is given, a custom component
    /// without predefined properties is created.
    void createComponents(BaseLib::ConfigTree const& config);

    /// The method creating phase properties based on config subtree.
    ///
    /// This method creates the properties of the Phase as defined in the
    /// prj-file. Only specified properties overwrite the default properties.
    void createProperties(BaseLib::ConfigTree const& config);

    /// The method initializing the defaule phase properties.
    /// Here, all for all properties listed in the Properties enumerator are
    /// initialized by mole average functions of value zero. However,
    /// 'special-default' properties are allowed to be set.
    void createDefaultProperties();

    /// A simple get-function for a component. The argument refers to the
    /// Index of the component in the components vector.
    Component& component(std::size_t const& index) const;

    /// A get-function for a component by component name.
    Component& component(std::string const& name) const;

    /// A get-function for a property. The argument refers to the name of the
    /// property.
    Property& property(PropertyEnum const& p) const;

    /// A get-function for retrieving the number of components in this phase.
    std::size_t numberOfComponents() const;
};

}  // namespace MaterialPropertyLib
