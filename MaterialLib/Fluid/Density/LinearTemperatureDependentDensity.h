/**
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file
 * \author: wenqing
 *
 * Created on August 10, 2016, 11:34 AM
 */

#pragma once

#include <string>

#include "BaseLib/ConfigTree.h"

#include "MaterialLib/Fluid/FluidProperty.h"

namespace MaterialLib
{
namespace Fluid
{
/// Linear temperature dependent density model.
class LinearTemperatureDependentDensity final : public FluidProperty
{
public:
    /**
     * @param rho0  \f$ \rho_0 \f$
     * @param T0    \f$ T_0 \f$
     * @param beta  \f$ \beta \f$
     */
    explicit LinearTemperatureDependentDensity(const double rho0, double T0,
                                               const double beta)
        : rho0_(rho0), temperature0_(T0), beta_(beta)
    {
    }

    /// Get model name.
    std::string getName() const override
    {
        return "Linear temperature dependent density";
    }

    /// Get density value.
    /// \param var_vals Variable values in an array. The order of its elements
    ///                 is given in enum class PropertyVariableType.
    double getValue(const ArrayType& var_vals) const override
    {
        const double T = var_vals[static_cast<int>(PropertyVariableType::T)];
        return rho0_ * (1 - beta_ * (T - temperature0_));
    }

    /// Get the partial differential of the density with respect to temperature.
    /// \param var_vals  Variable values  in an array. The order of its elements
    ///                   is given in enum class PropertyVariableType.
    /// \param var       Variable type.
    double getdValue(const ArrayType& var_vals,
                     const PropertyVariableType var) const override
    {
        (void)var_vals;
        if (var != PropertyVariableType::T)
        {
            return 0.0;
        }
        return -rho0_ * beta_;
    }

private:
    const double rho0_;          ///<  Reference density.
    const double temperature0_;  ///<  Reference temperature.
    const double beta_;          ///<  Parameter.
};

}  // namespace Fluid
}  // namespace MaterialLib
