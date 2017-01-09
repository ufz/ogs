/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file:   LinearTemperatureDependentDensity.h
 * \author: wenqing
 *
 * Created on August 10, 2016, 11:34 AM
 */

#ifndef LINEARTEMPERATUREDEPENDENTDENSITY_H
#define LINEARTEMPERATUREDEPENDENTDENSITY_H

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
        : _rho0(rho0), _temperature0(T0), _beta(beta)
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
        return _rho0 * (1 + _beta * (T - _temperature0));
    }

    /// Get the partial differential of the density with respect to temperature.
    /// \param var_vals  Variable values  in an array. The order of its elements
    ///                   is given in enum class PropertyVariableType.
    /// \param var       Variable type.
    double getdValue(const ArrayType& var_vals,
                     const PropertyVariableType var) const override
    {
        (void)var_vals;
        (void)var;
        return _rho0 * _beta;
    }

private:
    const double _rho0;          ///<  Reference density.
    const double _temperature0;  ///<  Reference temperature.
    const double _beta;          ///<  Parameter.
};

}  // end namespace
}  // end namespace

#endif /* LINEARTEMPERATUREDEPENDENTDENSITY_H */
