/**
 * @copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
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
/// (1 + \bar \alpha (C - C_{\text{ref}})) \f]
/// where
/// - \f$ \varrho_{\text{ref}}\f$ is the reference density
/// - \f$ \bar \alpha\f$ is the fluid density difference ratio
/// - \f$ C_{\text{ref}}\f$ is the reference concentration
class LinearConcentrationDependentDensity final : public FluidProperty
{
public:
    /**
     * @param reference_density  \f$\rho_0\f$
     * @param reference_concentration \f$C_0\f$
     * @param fluid_density_difference_ratio  \f$ \bar \alpha \f$ in reference
     * Coupled groundwater flow and transport: 2. Thermohaline and 3D convection
     * systems
     */
    explicit LinearConcentrationDependentDensity(
        const double reference_density, double reference_concentration,
        const double fluid_density_difference_ratio)
        : _reference_density(reference_density),
          _reference_concentration(reference_concentration),
          _fluid_density_difference_ratio(fluid_density_difference_ratio)
    {
    }

    /// Get model name.
    std::string getName() const override
    {
        return "Linear concentration dependent density";
    }

    /// Get density value.
    /// \param var_vals Variable values in an array. The order of its elements
    ///                 is given in enum class PropertyVariableType.
    double getValue(const ArrayType& var_vals) const override
    {
        const double C = var_vals[static_cast<int>(PropertyVariableType::C)];
        return _reference_density * (1 +
                                     _fluid_density_difference_ratio *
                                         (C - _reference_concentration));
    }

    /// Get the partial differential of the density with respect to
    /// concentration.
    /// \param var_vals  Variable values  in an array. The order of its elements
    ///                   is given in enum class PropertyVariableType.
    /// \param var       Variable type.
    double getdValue(const ArrayType& var_vals,
                     const PropertyVariableType var) const override
    {
        (void)var_vals;
        if (var != PropertyVariableType::C)
            return 0.0;
        return _reference_density * _fluid_density_difference_ratio;
    }

private:
    const double _reference_density;
    const double _reference_concentration;
    const double _fluid_density_difference_ratio;
};

}  // end namespace
}  // end namespace
