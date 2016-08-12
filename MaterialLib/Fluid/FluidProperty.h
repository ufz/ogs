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

namespace MaterialLib
{
namespace Fluid
{
/// Variable that determine the property.
enum class PropertyVariable
{
    T = 0,   ///< temperature
    pl = 1,  ///< pressure of the liquid phase (1st phase for some cases))
    pg = 2,  ///< pressure of the gas phase (2nd phase for some cases))
};

/// Base class of fluid density properties
class FluidProperty
{
public:
    FluidProperty() = default;

    /// Get model name.
    virtual std::string getName() const = 0;

    /// Get property value.
    /// \param var_vals Variable values in an array. The order of its elements
    ///                 is given in enum class PropertyVariable.
    virtual double getValue(const double /* var_vals*/[]) const = 0;

    /// Get the partial differential of the property value
    /// \param var_vals  Variable values  in an array. The order of its elements
    ///                   is given in enum class PropertyVariable.
    virtual double getdValue(const double /* var_vals*/[],
                             const PropertyVariable /* var */) const = 0;
};

}  // end namespace
}  // end namespace

#endif /* FLUIDPROPERTY_H */
