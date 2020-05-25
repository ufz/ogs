/**
 * @copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <string>

#include "BaseLib/ConfigTree.h"

#include "MaterialLib/Fluid/FluidProperty.h"

namespace MaterialLib
{
namespace Fluid
{
/// Linear concentration dependent density model.
/// \f[ \varrho = \varrho_{\text{ref}}
/// (1 + \bar \alpha (C - C_{\text{ref}}) + \bar \beta (p - p_{\text{ref}}) )
/// \f] where
/// - \f$ \varrho_{\text{ref}}\f$ is the reference density
/// - \f$ \bar \alpha\f$ is the fluid density concentration difference ratio
/// - \f$ C_{\text{ref}}\f$ is the reference concentration
/// - \f$ \bar \beta\f$ is the fluid density pressure difference ratio
/// - \f$ p_{\text{ref}}\f$ is the reference pressure
class LinearConcentrationAndPressureDependentDensity final
    : public FluidProperty
{
public:
    /**
     * @param reference_density  \f$\rho_0\f$
     * @param reference_concentration \f$C_0\f$
     * @param fluid_density_concentration_difference_ratio  \f$ \bar \alpha \f$
     * in reference
     * @param reference_pressure \f$p_0\f$
     * @param fluid_density_pressure_difference_ratio  \f$ \bar \beta \f$ in
     * reference Coupled groundwater flow and transport: 2. Thermohaline and 3D
     * convection systems
     */
    explicit LinearConcentrationAndPressureDependentDensity(
        const double reference_density, double reference_concentration,
        const double fluid_density_concentration_difference_ratio,
        double reference_pressure,
        const double fluid_density_pressure_difference_ratio)
        : reference_density_(reference_density),
          reference_concentration_(reference_concentration),
          fluid_density_concentration_difference_ratio_(
              fluid_density_concentration_difference_ratio),
          reference_pressure_(reference_pressure),
          fluid_density_pressure_difference_ratio_(
              fluid_density_pressure_difference_ratio)
    {
    }

    /// Get model name.
    std::string getName() const override
    {
        return "Linear concentration and pressure dependent density";
    }

    /// Get density value.
    /// \param var_vals Variable values in an array. The order of its elements
    ///                 is given in enum class PropertyVariableType.
    double getValue(const ArrayType& var_vals) const override
    {
        const double C = var_vals[static_cast<int>(PropertyVariableType::C)];
        const double p = var_vals[static_cast<int>(PropertyVariableType::p)];
        return reference_density_ *
               (1 +
                fluid_density_concentration_difference_ratio_ *
                    (C - reference_concentration_) +
                fluid_density_pressure_difference_ratio_ *
                    (p - reference_pressure_));
    }

    /// Get the partial differential of the density with respect to
    /// concentration or pressure.
    double getdValue(const ArrayType& /*var_vals*/,
                     const PropertyVariableType var) const override
    {
        if (var == PropertyVariableType::C)
        {
            return reference_density_ *
                   fluid_density_concentration_difference_ratio_;
        }
        if (var == PropertyVariableType::p)
        {
            return reference_density_ *
                   fluid_density_pressure_difference_ratio_;
        }
        return 0;
    }

private:
    const double reference_density_;
    const double reference_concentration_;
    const double fluid_density_concentration_difference_ratio_;
    const double reference_pressure_;
    const double fluid_density_pressure_difference_ratio_;
};

}  // namespace Fluid
}  // namespace MaterialLib
