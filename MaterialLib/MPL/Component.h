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
#pragma once

#include "BaseLib/Error.h"
#include "Property.h"

namespace MaterialPropertyLib
{
/// \brief This class defines components (substances).
/// \details The Component class is a base class used for not further specified
/// components. Components are specified by the property 'name'. For specified
/// components we derive special classes from this class (for clarity they are
/// located in the 'components' subfolder).
class Component
{
public:
    /// Default constructor of Component. This constructor is used
    /// when the component is not specified via the 'name'-tag.
    Component();

    /// Constructor for a custom component
    Component(std::string const& component_name,
              std::unique_ptr<PropertyArray>&& properties);
    virtual ~Component() = default;

    /// A get-function for retrieving a certain property.
    Property const& property(PropertyType const& /*p*/) const;

    Property const& operator[](PropertyType const& p) const;

    bool hasProperty(PropertyType const& p) const;

    template <typename T>
    T value(PropertyType const p) const
    {
        return property(p).template value<T>();
    }

    template <typename T>
    T value(PropertyType const p, VariableArray const& variable_array) const
    {
        return property(p).template value<T>(variable_array);
    }

    template <typename T>
    T dValue(PropertyType const p,
             VariableArray const& variable_array,
             Variable const variable) const
    {
        return property(p).template dValue<T>(variable_array, variable);
    }

    template <typename T>
    T d2Value(PropertyType const p,
              VariableArray const& variable_array,
              Variable const variable1,
              Variable const variable2) const
    {
        return property(p).template d2Value<T>(variable_array, variable1,
                                               variable2);
    }

    /// Short description of the component with its name.
    std::string description() const;

public:
    std::string const name;

protected:
    /// The property array of the component.
    PropertyArray properties_;
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
                                        bool& isCustomComponent);

template <typename Container>
void checkRequiredProperties(Component const& c,
                             Container const& required_properties)
{
    for (auto const& p : required_properties)
    {
        if (!c.hasProperty(p))
        {
            OGS_FATAL("The property '{:s}' is missing in the component '{:s}'.",
                      property_enum_to_string[p], c.name);
        }
    }
}

}  // namespace MaterialPropertyLib
