/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file
 *
 * Created on August 12, 2016, 3:34 PM
 */

#pragma once

#include <array>
#include <string>

#include "PropertyVariableType.h"

namespace MaterialLib
{
namespace Fluid
{
/// Base class of fluid properties
class FluidProperty
{
public:
    using ArrayType = std::array<double, PropertyVariableNumber>;

    virtual ~FluidProperty() = default;
    /// Get model name.
    virtual std::string getName() const = 0;

    /// Get property value.
    /// The argument is an array of variables. The order of its elements
    ///                 is given in enum class PropertyVariableType.
    virtual double getValue(const ArrayType& /* var_vals*/) const = 0;

    /// Get the partial differential of the property value
    /// The first argument is an array of variables, and the order of the array
    /// elements is given in enum class PropertyVariableType.
    /// The second argument is the variable type indicating which partial
    /// derivative to be calculated.
    virtual double getdValue(const ArrayType& /* var_vals*/,
                             const PropertyVariableType /* var */) const = 0;
};

}  // namespace Fluid
}  // namespace MaterialLib
