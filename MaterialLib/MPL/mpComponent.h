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

#include "BaseLib/ConfigTree.h"
#include "mpProperty.h"

namespace MaterialPropertyLib
{
/// \brief This class defines components (substances).
/// \details The Component class is a base class used for not further specified
/// components. Components are specified by the property 'name'. For specified
/// components we derive special classes from this class (for clarity they are
/// located in the 'components' subfolder).
class Component
{
protected:
    /// The property array of the component.
    PropertyArray _properties;

public:
    /// Default constructor of Component. This constructor is used
    /// when the component is not specified via the 'name'-tag.
    Component();

    /// Constructor for a custom component
    explicit Component(std::string const& component_name);
    virtual ~Component() = default;

    /// The method reads the 'properties' tag in the prj-file and creates
    /// component properties accordingly.
    ///
    /// First, a new property iy created based on the specified property type.
    /// Then, the property name is evaluated and the property is copied into the
    /// properties array.
    void createProperties(BaseLib::ConfigTree const& /*config*/);

    /// This method initializes the component property array with constant
    /// functions of value zero.
    void createDefaultProperties();
    /// A get-function for retrieving a cartain property.
    Property& property(PropertyEnum const& /*p*/) const;
};

/// Method for creating a new component based on the specified component name.
///
/// This function creates a new component based on the (optional) component name
/// that is given in the prj-file.
///
/// The method evaluates the string in the 'name'-object and calls the
/// constructors of the derived component classes (if found) or that of the base
/// class (if no name is specified).
std::unique_ptr<Component> newComponent(std::string const& component_name,
                                        bool& isNewComponent);

}  // namespace MaterialPropertyLib
