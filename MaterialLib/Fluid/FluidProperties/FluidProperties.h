/**
 *  \copyright
 *   Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *              Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file   FluidProperties.h
 *
 * Created on November 29, 2016, 2:23 PM
 */

#pragma once

#include <array>
#include <memory>

#include "MaterialLib/Fluid/FluidProperty.h"
#include "MaterialLib/Fluid/PropertyVariableType.h"

namespace MaterialLib
{
namespace Fluid
{
/// Fluid property type.
enum class FluidPropertyType
{
    Density = 0,
    Vicosity = 1,
    HeatCapacity = 2,
    ThermalConductivity = 3,
    number_of_property_types = 4  ///< Number of property types.
};

const unsigned FluidPropertyTypeNumber =
    static_cast<unsigned>(FluidPropertyType::number_of_property_types);

/// Base class of fluid properties.
class FluidProperties
{
public:
    typedef std::array<double, PropertyVariableNumber> ArrayType;

    FluidProperties(
        std::unique_ptr<MaterialLib::Fluid::FluidProperty>&& density,
        std::unique_ptr<MaterialLib::Fluid::FluidProperty>&& viscosity,
        std::unique_ptr<MaterialLib::Fluid::FluidProperty>&& heat_capacity,
        std::unique_ptr<MaterialLib::Fluid::FluidProperty>&&
            thermal_conductivity)
#ifndef MSVC
        // Up to the latest Visual Studio compiler, Visual Studio 2015
        // list initialization inside member initializer list or non-static data
        // member initializer is not implemented.
        // See: https://msdn.microsoft.com/en-us/library/dn793970.aspx
        : _property_models{{std::move(density), std::move(viscosity),
                            std::move(heat_capacity),
                            std::move(thermal_conductivity)}}
    {
    }
#else
    {
        // Move to the initialization list if MS visual studio compiler supports
        _property_models[static_cast<unsigned>(FluidPropertyType::Density)] =
            std::move(density);
        _property_models[static_cast<unsigned>(FluidPropertyType::Vicosity)] =
            std::move(viscosity);
        _property_models[static_cast<unsigned>(
            FluidPropertyType::HeatCapacity)] = std::move(heat_capacity);
        _property_models[static_cast<unsigned>(
            FluidPropertyType::ThermalConductivity)] =
            std::move(thermal_conductivity);
    }
#endif

    virtual ~FluidProperties() = default;

    /**
     *  Get the value of a Property.
     *  \param property_type   Property type.
     *  \param variable_values An array of the primary variables. The order of
     *                         its elements is temperature, pressure,
     *                         concentration, which is defined in enum class
     *                         PropertyVariableType.
     */
    virtual double getValue(const FluidPropertyType property_type,
                            const ArrayType& variable_values) const = 0;

    /**
     *  Get the partial differential of a property.
     *  \param property_type   Property type.
     *  \param variable_values An array of the primary variables. The order of
     *                         its elements is temperature, pressure,
     *                         concentration, which is defined in enum class
     *                         PropertyVariableType.
     *  \param variable_type   Variable type
     */
    virtual double getdValue(
        const FluidPropertyType property_type,
        const ArrayType& variable_values,
        const PropertyVariableType variable_type) const = 0;

protected:
#ifdef MSVC
    std::array<std::unique_ptr<FluidProperty>, FluidPropertyTypeNumber>
        _property_models;
#else
    /** Fluid property models.
     *  0: density;
     *  1: viscosity;
     *  2: specific heat capacity;
     *  3: thermal conductivity
     *
     *  The index is specified via enum class PropertyType.
     */
    const std::array<std::unique_ptr<FluidProperty>, FluidPropertyTypeNumber>
        _property_models;
#endif
};

}  // end namespace
}  // end namespace
