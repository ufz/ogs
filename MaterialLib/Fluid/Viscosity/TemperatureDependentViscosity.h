/**
 *  \brief A linear temperature dependent viscosity model.
 *
 *  \copyright
 *   Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *              Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file:   TemperatureDependentViscosity.h
 *
 *  Created on August 10, 2016, 10:45 AM
 */

#ifndef TEMPERATUREDEPENDENTVISCOSITY_H
#define TEMPERATUREDEPENDENTVISCOSITY_H

#include <string>
#include <cmath>

#include "MaterialLib/Fluid/FluidProperty.h"

namespace MaterialLib
{
namespace Fluid
{
/**
 *  Class of a temperature dependent viscosity model, which provides a good fit
 *  over temperature variations of 150 °C.
 *
 *   \f[
 *        \mu = \mu_0\, exp(-\frac{T-T_c}{T_v})
 *   \f]
 *   where
 *    \f{eqnarray*}{
 *       &\mu_0:&  \mbox{reference viscosity at the temperature of} T_c\\
 *       &T:&      \mbox{temperature,}\\
 *       &T_c:&    \mbox{reference temperature,}\\
 *       &T_v:&    \mbox{constant}\\
 *    \f}
 */
class TemperatureDependentViscosity final : public FluidProperty
{
public:
    /**
     *  @param mu0   \f$ \mu_0 \f$
     *  @param T_c   \f$ T_c \f$
     *  @param T_v   \f$ T_v \f$
     */
    explicit TemperatureDependentViscosity(const double mu0,
                                           const double T_c,
                                           const double T_v)
        : _mu0(mu0), _temperature_c(T_c), _temperature_v(T_v)
    {
    }

    /// Get model name.
    std::string getName() const override
    {
        return "Temperature dependent viscosity";
    }

    /** Get viscosity value.
     * \param var_vals Variable values in an array. The order of its elements
     *                 is given in enum class PropertyVariableType.
     */
    double getValue(const ArrayType& var_vals) const override
    {
        const double T = var_vals[static_cast<int>(PropertyVariableType::T)];
        return _mu0 * std::exp(-(T - _temperature_c) / _temperature_v);
    }

    /**
     *  Get the partial differential of the viscosity with respect to
     * temperature.
     * \param var_vals  Variable values  in an array. The order of its elements
     *                  is given in enum class PropertyVariableType.
     * \param var       Variable type.
     */
    double getdValue(const ArrayType& var_vals,
                     const PropertyVariableType var) const override
    {
        (void)var;
        const double T = var_vals[static_cast<int>(PropertyVariableType::T)];
        return -_mu0 * std::exp(-(T - _temperature_c) / _temperature_v);
    }

private:
    const double _mu0;            ///<  Inital viscosity.
    const double _temperature_c;  ///<  Reference temperature 1.
    const double _temperature_v;  ///<  Reference temperature 2.
};

}  // end namespace
}  // end namespace

#endif /* TEMPERATUREDEPENDENTVISCOSITY_H */
