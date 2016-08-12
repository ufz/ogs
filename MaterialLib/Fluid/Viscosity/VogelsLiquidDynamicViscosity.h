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

#include "ViscosityType.h"

namespace MaterialLib
{
namespace Fluid
{
class VogelsLiquidDynamicViscosity
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
        assert(mat_id < 3);  // So far only three sets of data are given.
    }

    /// Get model name.
    std::string getName() const
    {
        return "Liquid viscosity by Vogel's equation";
    }

    ViscosityType getType() const { return ViscosityType::VOGEL; }
    /// Get viscosity value
    /// \param T Temperature
    double getValue(const double T) const
    {
        return 1.e-3 * std::exp(_abc[0] + _abc[1] / (_abc[2] + T));
    }

    /// Get the derivative of viscosity
    /// \param T Temperature
    double getdValue(const double T) const
    {
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

    /// Address of the data array in _constants for given liquid type.
    const std::vector<double>& _abc;
};

}  // end namespace
}  // end namespace
#endif
