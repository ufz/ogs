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
/// Constant fluid properties
class ConstantFluidProperty final : public FluidProperty
{
public:
    explicit ConstantFluidProperty(const double value) : FluidProperty(),_value(value)
    {
    }

    /// Get model name.
    std::string getName() const override
    {
        return "Constant";
    }

    /// Get property value.
    /// \param var_vals Variable values in an array. The order of its elements
    ///                 is given in enum class PropertyVariable.
    double getValue(const double var_vals[]) const override
    {
        (void) var_vals;
        return _value;
    }

    /// Get the partial differential of the property value
    /// \param var_vals  Variable values  in an array. The order of its elements
    ///                  is given in enum class PropertyVariable.
    /// \param var       Variable type.
    double getdValue(const double var_vals[],
            const PropertyVariable var) const override
    {
        (void) var_vals;
        (void) var;
        return 0.;
    }

private:
    const double _value;
};

}  // end namespace
}  // end namespace
#endif /* CONSTANTFLUIDPROPERTY_H */
