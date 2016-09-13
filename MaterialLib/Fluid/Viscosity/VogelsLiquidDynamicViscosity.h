/*!
   \file  VogelsLiquidDynamicViscosity.h
   \brief Declaration of class for the pressure dependent viscosity model.

   \author Wenqing Wang
   \date Jan 2015

   \copyright
    Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
               Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
 */
#ifndef VOGELS_LIQUID_DYNAMIC_VISCOSITY_H_
#define VOGELS_LIQUID_DYNAMIC_VISCOSITY_H_

#include <string>
#include <vector>
#include <cmath>

#include "MaterialLib/Fluid/FluidProperty.h"

namespace MaterialLib
{
namespace Fluid
{
/// A temperature dependent viscosity model.
class VogelsLiquidDynamicViscosity : public FluidProperty
{
public:
    /*!
         \brief Viscosity defined by \f$e^{A+\dfrac{B}{C+T}}\f$
         \param mat_id  0: Water
                        1: Carbon dioxide
                        2: Methane
     */
    VogelsLiquidDynamicViscosity(const unsigned mat_id)
        : _abc(_constants[mat_id])
    {
    }

    /// Get model name.
    virtual std::string getName() const final
    {
        return "Liquid viscosity by Vogel's equation";
    }

    /// Get viscosity value.
    /// \param var_vals Variable values in an array. The order of its elements
    ///                 is given in enum class PropertyVariable.
    virtual double getValue(const double var_vals[]) const final
    {
        return 1.e-3 *
               std::exp(_abc[0] +
                        _abc[1] /
                            (_abc[2] +
                             var_vals[static_cast<int>(PropertyVariable::T)]));
    }

    /// Get the partial differential of the viscosity with respect to
    /// temperature.
    /// \param var_vals  Variable values  in an array. The order of its elements
    ///                   is given in enum class PropertyVariable.
    virtual double getdValue(const double var_vals[],
                             const PropertyVariable /* var */) const final
    {
        const double T = var_vals[static_cast<int>(PropertyVariable::T)];
        const double f_buff = _abc[1] / (_abc[2] + T);
        return -1.e-3 * f_buff * std::exp(_abc[0] + f_buff) / (_abc[2] + T);
    }

private:
    ///  Parameters A, B, C.
    /// <a href="ddbst">
    /// http://ddbonline.ddbst.de/VogelCalculation/VogelCalculationCGI.exe</a>
    const std::vector<std::vector<double>> _constants = {
        {-3.7188, 578.919, -137.546},  // Water
        {-24.0592, 28535.2, 1037.41},  // Carbon dioxide
        {-25.5947, 25392, 969.306}};   // Methane

    /// Parameter vector in _constants for a given liquid type.
    const std::vector<double>& _abc;
};

}  // end namespace
}  // end namespace
#endif
