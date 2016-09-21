/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file:   FluidProperty.h
 *
 * Created on August 12, 2016, 3:34 PM
 */

#ifndef FLUIDPROPERTY_H
#define FLUIDPROPERTY_H

#include <string>
#include <array>

namespace MaterialLib
{
namespace Fluid
{
/// Variable that determine the property.
enum class PropertyVariableType
{
    T = 0,///< temperature.
    pl = 1,///< pressure of the liquid phase (1st phase for some cases).
    pg = 2,///< pressure of the gas phase (2nd phase for some cases).
    number_of_variables = 3 ///< Number of property variables.
};

const unsigned PropertyVariableNumber
        = static_cast<unsigned> (PropertyVariableType::number_of_variables);

/// Base class of fluid density properties
class FluidProperty
{
public:

    typedef std::array<double,PropertyVariableNumber> ArrayType;

    virtual ~FluidProperty()
    {
    }

    /// Get model name.
    virtual std::string getName() const = 0;

    /// Get property value.
    /// The argument is an array of variables. The order of its elements
    ///                 is given in enum class PropertyVariableType.
    virtual double getValue(const ArrayType& /* var_vals*/) const = 0;

    /// Get the partial differential of the property value
    /// The first argument is an array of variables, and the order of the array
    /// elements is given in enum class PropertyVariableType.
    /// The second argument is the variable type indicating which partial derivative
    /// to be calculated.
    virtual double getdValue(const ArrayType& /* var_vals*/,
            const PropertyVariableType /* var */) const = 0;
};

}  // end namespace
}  // end namespace

#endif /* FLUIDPROPERTY_H */
