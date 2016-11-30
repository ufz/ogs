/**
 *  \copyright
 *   Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *              Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file   TemperaturePressureConcentrationDependentFluidProperties.h
 *
 * Created on November 29, 2016, 3:19 PM
 */

#ifndef OGS_TEMPERATURE_PRESSURE_CONCENTRATION_DEPENDENT_FLUID_PROPERTIES_H
#define OGS_TEMPERATURE_PRESSURE_CONCENTRATION_DEPENDENT_FLUID_PROPERTIES_H

#include "FluidProperties.h"
#include "MaterialLib/Fluid/FluidProperty.h"

namespace MaterialLib
{
namespace Fluid
{
class FluidProperty;

///  A class contains density, viscosity, heat_capacity and thermal_conductivity
///  models, which are all functions of temperature, pressure and concentration.
class TemperaturePressureConcentrationDependentFluidProperties final
    : public FluidProperties
{
public:
    TemperaturePressureConcentrationDependentFluidProperties(
        std::unique_ptr<MaterialLib::Fluid::FluidProperty>&& density,
        std::unique_ptr<MaterialLib::Fluid::FluidProperty>&& viscosity,
        std::unique_ptr<MaterialLib::Fluid::FluidProperty>&& heat_capacity,
        std::unique_ptr<MaterialLib::Fluid::FluidProperty>&&
            thermal_conductivity)
        : FluidProperties(std::move(density), std::move(viscosity),
                          std::move(heat_capacity),
                          std::move(thermal_conductivity))
    {
    }

    /**
     *  Get the value of a Property.
     *  \param property_type   Property type.
     *  \param variable_values An array of variables. The order of its elements
     *                         is temperature, pressure, concentration, which is
     *                         defined in enum class PropertyVariableType.
     */
    double getValue(const FluidPropertyType property_type,
                    const ArrayType& variable_values) const override
    {
        return _property_models[static_cast<unsigned>(property_type)]->getValue(
            variable_values);
    }

    /**
     *  Get the partial differential of a property.
     *  \param property_type   Property type.
     *  \param variable_values An array of variables. The order of its elements
     *                         is temperature, pressure, concentration, which is
     *                         defined in enum class PropertyVariableType.
     *  \param variable_type   Variable type
     */
    double getdValue(const FluidPropertyType property_type,
                     const ArrayType& variable_values,
                     const PropertyVariableType variable_type) const override
    {
        return _property_models[static_cast<unsigned>(property_type)]
            ->getdValue(variable_values, variable_type);
    }
};

}  // end namespace
}  // end namespace
#endif /* OGS_TEMPERATURE_PRESSURE_CONCENTRATION_DEPENDENT_FLUID_PROPERTIES_H \
          */
