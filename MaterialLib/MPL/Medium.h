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
 *
 */
#pragma once

#include <memory>
#include <string>
#include <vector>

#include "Phase.h"

namespace MaterialPropertyLib
{
class Property;
}
namespace MaterialPropertyLib
{
/// This class is for material objects on the Medium scale.
///
/// A Medium consists of an arbitrarily long vector of phases and an array of
/// properties.
class Medium final
{
public:
    Medium(std::vector<std::unique_ptr<Phase>>&& phases,
           std::unique_ptr<PropertyArray>&& properties);

    /// A get-function for a particular phase. The ul argument specifies the
    /// index within the _phases vector.
    Phase const& phase(std::size_t index) const;
    /// A get-function for a particular phase by phase name.
    Phase const& phase(std::string const& phase_name) const;
    /// A get-function for a property. The argument refers to the name of the
    /// property.
    Property const& property(PropertyType const& p) const;

    /// A simple get-function for retrieving the number of phases the medium
    /// consists of.
    std::size_t numberOfPhases() const;

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
             Variables const variables) const
    {
        return property(p).template dValue<T>(variable_array, variables);
    }

    template <typename T>
    T d2Value(PropertyType const p,
              VariableArray const& variable_array,
              Variables const variables1,
              Variables const variables2) const
    {
        return property(p).template d2Value<T>(variable_array, variables1,
                                               variables2);
    }

private:
    /// The vector that holds the phases.
    std::vector<std::unique_ptr<Phase>> const _phases;
    /// The array that holds the medium properties.
    ///
    /// Currently, these defaults is the volume fraction average.
    ///
    /// Most properties are fine with the volume fraction average, but
    /// special-defaults are allowed as well...
    PropertyArray _properties;
};
}  // namespace MaterialPropertyLib
