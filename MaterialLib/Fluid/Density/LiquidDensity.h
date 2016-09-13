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
class LiquidDensity : public FluidProperty
{
public:
    /// \param config  ConfigTree object which contains the input data
    ///                 including  <type>fluid</type> and it has
    ///                 a tag of <density>
    LiquidDensity(BaseLib::ConfigTree const& config)
        :  //! \ogs_file_param{material__fluid__density__liquid_density__beta}
    _beta(config.getConfigParameter<double>("beta")),
          //! \ogs_file_param{material__fluid__density__liquid_density__rho0}
    _rho0(config.getConfigParameter<double>("rho0")),
          //! \ogs_file_param{material__fluid__density__liquid_density__temperature0}
    _temperature0(config.getConfigParameter<double>("temperature0")),
          //! \ogs_file_param{material__fluid__density__liquid_density__p0}
    _p0(config.getConfigParameter<double>("p0")),
          //! \ogs_file_param{material__fluid__density__liquid_density__bulk_modulus}
    _bulk_moudlus(config.getConfigParameter<double>("bulk_modulus"))
    {
    }

    /// Get density model name.
    virtual std::string getName() const final
    {
        return "Temperature or pressure dependent liquid density";
    }

    /// Get density value.
    /// \param var_vals Variable values in an array. The order of its elements
    ///                 is given in enum class PropertyVariable.
    virtual double getValue(const double var_vals[]) const final
    {
        const double T = var_vals[static_cast<int>(PropertyVariable::T)];
        const double p = var_vals[static_cast<int>(PropertyVariable::pl)];
        return _rho0 / (1. + _beta * (T - _temperature0)) /
               (1. - (p - _p0) / _bulk_moudlus);
    }

    /// Get the partial differential of the density with respect to temperature
    /// or liquid pressure.
    /// \param var_vals  Variable values  in an array. The order of its elements
    ///                   is given in enum class PropertyVariable.
    virtual double getdValue(const double var_vals[],
                             const PropertyVariable var) const final
    {
        assert(var == PropertyVariable::T || var == PropertyVariable::pl);

        const int func_id = (var == PropertyVariable::T) ? 0 : 1;
        const double T = var_vals[static_cast<int>(PropertyVariable::T)];
        const double p = var_vals[static_cast<int>(PropertyVariable::pl)];
        return (this->*_derivative_functions[func_id])(T, p);
    }

private:
    /// Volumetric temperature expansion coefficient.
    double _beta;

    /// Initial density.
    double _rho0;

    /// Initial temperature.
    double _temperature0;

    /// Initial pressure.
    double _p0;

    /// Bulk modulus.
    double _bulk_moudlus;

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
               (1. - (p - _p0) / _bulk_moudlus);
    }

    /**
     *  Calculate the derivative of density of fluids with respect to pressure.
     *   \param T    Temperature (K).
     *   \param p    Pressure (Pa).
     */
    double dLiquidDensity_dp(const double T, const double p) const
    {
        const double fac_p = 1. - (p - _p0) / _bulk_moudlus;
        return _rho0 / (1. + _beta * (T - _temperature0)) /
               (fac_p * fac_p * _bulk_moudlus);
    }

    typedef double (LiquidDensity::*DerivativeFunctionPointer)(const double,
                                                       const double) const;

    /// An array of pointer to derivative functions.
    static DerivativeFunctionPointer _derivative_functions[2];
};

LiquidDensity::DerivativeFunctionPointer LiquidDensity::_derivative_functions[2]
        = {&LiquidDensity::dLiquidDensity_dT,&LiquidDensity::dLiquidDensity_dp};

}  // end of namespace
}  // end of namespace

#endif /* LIQUIDDENSITY_H */
