/**
 *  \copyright
 *   Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *              Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file   FluidPropertiesWithDensityDependentModels.cpp
 *
 * Created on November 29, 2016, 3:19 PM
 */

#include "FluidPropertiesWithDensityDependentModels.h"

#include <string>

#include "MaterialLib/Fluid/FluidProperty.h"
#include "FluidProperties.h"

namespace MaterialLib
{
namespace Fluid
{
FluidPropertiesWithDensityDependentModels::
    FluidPropertiesWithDensityDependentModels(
        std::unique_ptr<MaterialLib::Fluid::FluidProperty>&& density,
        std::unique_ptr<MaterialLib::Fluid::FluidProperty>&& viscosity,
        std::unique_ptr<MaterialLib::Fluid::FluidProperty>&& heat_capacity,
        std::unique_ptr<MaterialLib::Fluid::FluidProperty>&&
            thermal_conductivity,
        const bool is_viscosity_density_dependent,
        const bool is_heat_capacity_dependent,
        const bool is_thermal_conductivity)
    : FluidProperties(std::move(density), std::move(viscosity),
                      std::move(heat_capacity),
                      std::move(thermal_conductivity)),
      _is_density_dependent{{false, is_viscosity_density_dependent,
                             is_heat_capacity_dependent,
                             is_thermal_conductivity}}
{
}

double FluidPropertiesWithDensityDependentModels::getValue(
    const FluidPropertyType property_type,
    const ArrayType& variable_values) const
{
    ArrayType var_vals = variable_values;
    if (_is_density_dependent[static_cast<unsigned>(property_type)])
    {
        var_vals[static_cast<unsigned>(PropertyVariableType::rho)] =
            _property_models[static_cast<unsigned>(FluidPropertyType::Density)]
                ->getValue(variable_values);
    }
    return _property_models[static_cast<unsigned>(property_type)]->getValue(
        var_vals);
}

double FluidPropertiesWithDensityDependentModels::getdValue(
    const FluidPropertyType property_type,
    const ArrayType& variable_values,
    const PropertyVariableType variable_type) const
{
    if (_is_density_dependent[static_cast<unsigned>(property_type)])
    {
        if (variable_type == PropertyVariableType::T)
        {
            return compute_df_drho_drho_dT(property_type, variable_values);
        }
        if (variable_type == PropertyVariableType::p)
        {
            return compute_df_drho_drho_dp(property_type, variable_values);
        }
    }
    else
    {
        return _property_models[static_cast<unsigned>(property_type)]
            ->getdValue(variable_values, variable_type);
    }

    return 0.0;
}

double FluidPropertiesWithDensityDependentModels::compute_df_drho_drho_dT(
    const FluidPropertyType property_type,
    const ArrayType& variable_values) const
{
    const auto& fluid_density_model =
        _property_models[static_cast<unsigned>(FluidPropertyType::Density)];
    const double drho_dT = fluid_density_model->getdValue(
        variable_values, PropertyVariableType::T);

    const double density_value = fluid_density_model->getValue(variable_values);

    ArrayType var_vals = variable_values;
    var_vals[static_cast<unsigned>(PropertyVariableType::rho)] = density_value;
    const auto& fluid_property_model =
        _property_models[static_cast<unsigned>(property_type)];

    // return d()/dT + d ()/drho * drho/dT
    return fluid_property_model->getdValue(var_vals, PropertyVariableType::T) +
           fluid_property_model->getdValue(var_vals,
                                           PropertyVariableType::rho) *
               drho_dT;
}

double FluidPropertiesWithDensityDependentModels::compute_df_drho_drho_dp(
    const FluidPropertyType property_type,
    const ArrayType& variable_values) const
{
    const auto& fluid_density_model =
        _property_models[static_cast<unsigned>(FluidPropertyType::Density)];

    const double drho_dp = fluid_density_model->getdValue(
        variable_values, PropertyVariableType::p);

    const double density_value = fluid_density_model->getValue(variable_values);
    ArrayType var_vals = variable_values;
    var_vals[static_cast<unsigned>(PropertyVariableType::rho)] = density_value;

    // return  d ()/drho * drho/dp
    return _property_models[static_cast<unsigned>(property_type)]->getdValue(
               var_vals, PropertyVariableType::rho) *
           drho_dp;
}

}  // end namespace
}  // end namespace
