/**
 *  \copyright
 *   Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *              Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file   CompositeFluidProperty.h
 *
 * Created on November 29, 2016, 2:23 PM
 */

#ifndef OGS_COMPOSITE_FLUID_PROPERTY_H
#define OGS_COMPOSITE_FLUID_PROPERTY_H

#include <array>

#include "PropertyVariableType.h"

namespace MaterialLib
{
namespace Fluid
{
/// Property type.
enum class PropertyType
{
    Density,
    Vicosity,
    HeatCapacity,
    ThermalConductivity
};

class CompositeFluidProperty
{
public:
    typedef std::array<double, PropertyVariableNumber> ArrayType;

    virtual ~CompositeFluidProperty(){};

    /**
     *  Get the value of a Property.
     *  \param property_type   Property type.
     *  \param variable_values An array of variables. The order of its elements
     *                         is given in enum class PropertyVariableType.
     */
    virtual double getValue(const PropertyType property_type,
                            const ArrayType& variable_values) const = 0;

    /**
     *  Get the partial differential of a property.
     *  \param property_type   Property type.
     *  \param variable_values An array of variables. The order of its elements
     *                         is given in enum class PropertyVariableType.
     *  \param variable_type   Variable type
     */
    virtual double getdValue(
        const PropertyType property_type,
        const ArrayType& variable_values,
        const PropertyVariableType variable_type) const = 0;
};

}  // end namespace
}  // end namespace
#endif /* OGS_COMPOSITE_FLUID_PROPERTY_H */
