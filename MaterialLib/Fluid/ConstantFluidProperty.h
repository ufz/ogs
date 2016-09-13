/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file:   ConstantFluidProperty.h
 *
 * Created on August 15, 2016, 12:11 PM
 */

#ifndef CONSTANTFLUIDPROPERTY_H
#define CONSTANTFLUIDPROPERTY_H

#include "FluidProperty.h"

namespace MaterialLib
{
namespace Fluid
{
/// Constant fluid density properties
class ConstantFluidProperty : public FluidProperty
{
public:
    ConstantFluidProperty(const double value) : FluidProperty(), _value(value)
    {
    }

    /// Get model name.
    virtual std::string getName() const final { return "constant"; }
    /// Get property value.
    /// \param var_vals Variable values in an array. The order of its elements
    ///                 is given in enum class PropertyVariable.
    virtual double getValue(const double /* var_vals*/[]) const final
    {
        return _value;
    }

    /// Get the partial differential of the property value
    /// \param var_vals  Variable values  in an array. The order of its elements
    ///                   is given in enum class PropertyVariable.
    virtual double getdValue(const double /* var_vals*/[],
                             const PropertyVariable /* var */) const final
    {
        return 0.;
    }

private:
    double _value;
};

}  // end namespace
}  // end namespace
#endif /* CONSTANTFLUIDPROPERTY_H */
