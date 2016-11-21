/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \author wenqing
 * \file   LiquidDensity.h
 *
 * Created on August 4, 2016, 10:15 AM
 */

#ifndef LIQUIDDENSITY_H
#define LIQUIDDENSITY_H

#include <vector>

#include "BaseLib/Error.h"

#include "BaseLib/ConfigTree.h"

#include "MaterialLib/Fluid/FluidProperty.h"

namespace MaterialLib
{
namespace Fluid
{
/**
 *  Class of a generic density model for liquid fluids varying with temperature
 *  or pressure.
 *
 *  The formula is given on
 *  <a href="engineeringtoolbox">
 * http://www.engineeringtoolbox.com/fluid-density-temperature-pressure-d_309.html</a>
 *   which reads
 *   \f[
 *        \rho_l = \rho_0/(1+\beta(T-T_0))/(1-(p-p_0)/E)
 *   \f]
 *   where
 *    \f{eqnarray*}{
 *       &\rho_l:& \mbox{liquid density,}\\
 *       &\rho_0:& \mbox{initial liquid density,}\\
 *       &\beta: & \mbox{volumetric temperature expansion coefficient,}\\
 *       &T:&      \mbox{temperature,}\\
 *       &T_0:&    \mbox{initial temperature,}\\
 *       &p:&      \mbox{pressure,}\\
 *       &p_0:&    \mbox{initial pressure,}\\
 *       &E:&      \mbox{bulk modulus fluid elasticity.}\\
 *    \f}
 */
class LiquidDensity final : public FluidProperty
{
public:
    /**
     * @param beta \$f \beta \f$
     * @param rho  \$f \rho_0 \f$
     * @param T0   \$f T_0 \f$
     * @param p0   \$f p_0 \f$
     * @param E    \$f E \f$
     */
    explicit LiquidDensity(const double beta, const double rho0,
                           const double T0, const double p0, const double E)
        : _beta(beta), _rho0(rho0), _temperature0(T0), _p0(p0), _bulk_modulus(E)
    {
    }

    /// Get density model name.
    std::string getName() const override
    {
        return "Temperature or pressure dependent liquid density";
    }

    /// Get density value.
    /// \param var_vals Variable values in an array. The order of its elements
    ///                 is given in enum class PropertyVariableType.
    double getValue(const ArrayType& var_vals) const override
    {
        const double T = var_vals[static_cast<int>(PropertyVariableType::T)];
        const double p = var_vals[static_cast<int>(PropertyVariableType::pl)];
        return _rho0 / (1. + _beta * (T - _temperature0)) /
               (1. - (p - _p0) / _bulk_modulus);
    }

    /// Get the partial differential of the density with respect to temperature
    /// or liquid pressure.
    /// \param var_vals  Variable values  in an array. The order of its elements
    ///                   is given in enum class PropertyVariableType.
    /// \param var       Variable type.
    double getdValue(const ArrayType& var_vals,
                     const PropertyVariableType var) const override
    {
        const double T = var_vals[static_cast<int>(PropertyVariableType::T)];
        const double p = var_vals[static_cast<int>(PropertyVariableType::pl)];
        switch (var)
        {
            case PropertyVariableType::T:
                return dLiquidDensity_dT(T, p);
            case PropertyVariableType::pl:
                return dLiquidDensity_dp(T, p);
            default:
                return 0.;
        }
    }

private:
    /// Volumetric temperature expansion coefficient.
    const double _beta;

    /// Initial density.
    const double _rho0;

    /// Initial temperature.
    const double _temperature0;

    /// Initial pressure.
    const double _p0;

    /// Bulk modulus.
    const double _bulk_modulus;

    /**
     *  Calculate the derivative of density of fluids with respect to
     * temperature.
     *   \param T    Temperature (K).
     *   \param p    Pressure (Pa).
     */
    double dLiquidDensity_dT(const double T, const double p) const
    {
        const double fac_T = 1. + _beta * (T - _temperature0);
        return -_beta * _rho0 / (fac_T * fac_T) /
               (1. - (p - _p0) / _bulk_modulus);
    }

    /**
     *  Calculate the derivative of density of fluids with respect to pressure.
     *   \param T    Temperature (K).
     *   \param p    Pressure (Pa).
     */
    double dLiquidDensity_dp(const double T, const double p) const
    {
        const double fac_p = 1. - (p - _p0) / _bulk_modulus;
        return _rho0 / (1. + _beta * (T - _temperature0)) /
               (fac_p * fac_p * _bulk_modulus);
    }
};

}  // end of namespace
}  // end of namespace

#endif /* LIQUIDDENSITY_H */
