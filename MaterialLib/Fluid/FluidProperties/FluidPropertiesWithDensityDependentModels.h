/**
 *  \copyright
 *   Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *              Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file   FluidPropertiesWithDensityDependentModels.h
 *
 * Created on November 29, 2016, 4:18 PM
 */

#ifndef OGS_FLUID_PROPERTIES_WITH_DENSITY_DEPENDENT_MODELS_H
#define OGS_FLUID_PROPERTIES_WITH_DENSITY_DEPENDENT_MODELS_H

#include <atomic>

#include "FluidProperties.h"

namespace MaterialLib
{
namespace Fluid
{
class FluidProperty;

/**
 * A class contains density, viscosity, heat_capacity and thermal_conductivity
 *   models, which are all functions of temperature, density and concentration.
 *
 *   \attention The density must be computed firstly, i.e. first call
 *             getValue(PropertyType::density, variable_values).
 */
class FluidPropertiesWithDensityDependentModels final : public FluidProperties
{
public:
    FluidPropertiesWithDensityDependentModels(
        std::unique_ptr<MaterialLib::Fluid::FluidProperty>&& density,
        std::unique_ptr<MaterialLib::Fluid::FluidProperty>&& viscosity,
        std::unique_ptr<MaterialLib::Fluid::FluidProperty>&& heat_capacity,
        std::unique_ptr<MaterialLib::Fluid::FluidProperty>&&
            thermal_conductivity,
        const bool is_viscosity_density_dependent,
        const bool is_heat_capacity_dependent,
        const bool is_thermal_conductivity);

    /**
     *  Get the value of a Property.
     *  \param property_type   Property type.
     *  \param variable_values An array of the primary variables. The order of
     *                         its elements is temperature, pressure,
     *                         concentration, which is defined in enum class
     *                         PropertyVariableType.
     *  \attention The density must be computed firstly, i.e. first call
     *             getValue(PropertyType::density, variable_values).
     */
    double getValue(const FluidPropertyType property_type,
                    const ArrayType& variable_values) const override;

    /**
     *  Get the partial differential of a property.
     *  \param property_type   Property type.
     *  \param variable_values An array of the primary variables. The order of
     *                         its elements is temperature, pressure,
     *                         concentration, which is defined in enum class
     *                         PropertyVariableType.
     *  \param variable_type   Variable type
     */
    double getdValue(const FluidPropertyType property_type,
                     const ArrayType& variable_values,
                     const PropertyVariableType variable_type) const override;

private:
    /// Compute df/dT for f(T, rho) with rho(T, p)
    double compute_df_drho_drho_dT(const double density_value,
                                   const FluidPropertyType property_type,
                                   const ArrayType& variable_values) const;
    /// Compute df/dp for f(T, rho) with rho(T, p)
    double compute_df_drho_drho_dp(const double density_value,
                                   const FluidPropertyType property_type,
                                   const ArrayType& variable_values) const;

    std::array<bool, FluidPropertyTypeNumber> _is_density_depedent;
};

}  // end namespace
}  // end namespace
#endif /* OGS_FLUID_PROPERTIES_WITH_DENSITY_DEPENDENT_MODELS_H */
