/**
 *  \brief Declaration of class for the pressure dependent viscosity model.
 *
 *  \copyright
 *   Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *              Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *   \file  LinearPressureDependentViscosity.h
 */
#ifndef LINEAR_PRESSURE_DEPENDENT_VISCOSITY_H_
#define LINEAR_PRESSURE_DEPENDENT_VISCOSITY_H_

#include <string>

#include "BaseLib/ConfigTree.h"

#include "MaterialLib/Fluid/FluidProperty.h"

namespace MaterialLib
{
namespace Fluid
{
/**
 *  Class of a linear pressure dependent viscosity model.
 *   \f[
 *        \mu = \mu_0\,(1 + \gamma (p -p_0) )
 *   where
 *    \f{eqnarray*}{
 *       &\mu_0:&  \mbox{reference viscosity,}\\
 *       &\gamma:& \mbox{parameter,}\\
 *       &p:&      \mbox{pressure,}\\
 *       &p_0:&    \mbox{initial pressure.}\\
 *    \f}
 */
class LinearPressureDependentViscosity final : public FluidProperty
{
public:
    /** \param parameters An array contains the three parameters:
     *                     [0] $f\mu_0$f
     *                     [1] $fp_0$f
     *                     [2] $f\gamma$f
     */
    explicit LinearPressureDependentViscosity(
        std::array<double, 3> const& parameters)
        : _mu0(parameters[0]), _p0(parameters[1]), _gamma(parameters[2])
    {
    }

    /// Get model name.
    std::string getName() const override
    {
        return "Linear pressure dependent viscosity";
    }

    /** Get viscosity value.
     * \param var_vals Variable values in an array. The order of its elements
     *                 is given in enum class PropertyVariableType.
     */
    double getValue(const ArrayType& var_vals) const override
    {
        const double p = var_vals[static_cast<int>(PropertyVariableType::pl)];
        return _mu0 * (1 + _gamma * (std::max(p, 0.) - _p0));
    }

    /** Get the partial differential of the viscosity with respect to pressure.
     * \param var_vals  Variable values  in an array. The order of its elements
     *                  is given in enum class PropertyVariableType.
     * \param var       Variable type.
     */
    double getdValue(const ArrayType& var_vals,
                     const PropertyVariableType var) const override
    {
        (void)var_vals;
        (void)var;
        return _mu0 * _gamma;
    }

private:
    const double _mu0;    ///<  Reference viscosity.
    const double _p0;     ///<  Reference pressure.
    const double _gamma;  ///<  Parameter.
};

}  // end namespace
}  // end namespace
#endif
