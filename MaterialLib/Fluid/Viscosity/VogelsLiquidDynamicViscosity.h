/**
 *  \brief Declaration of class for the pressure dependent viscosity model.
 *
 *  \copyright
 *   Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *              Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  \file  VogelsLiquidDynamicViscosity.h
 *
 */

#ifndef VOGELS_LIQUID_DYNAMIC_VISCOSITY_H_
#define VOGELS_LIQUID_DYNAMIC_VISCOSITY_H_

#include <string>
#include <vector>
#include <array>
#include <cmath>

#include "MaterialLib/Fluid/FluidProperty.h"

namespace MaterialLib
{
namespace Fluid
{
/** A temperature dependent viscosity model.
 * <a href="ddbst"> http://ddbonline.ddbst.de/VogelCalculation/VogelCalculationCGI.exe</a>
 */
template <typename VogelsConstants>
class VogelsLiquidDynamicViscosity final : public FluidProperty
{
public:
    /**
     *  \brief Viscosity defined by \f$10^3 \, e^{A+\frac{B}{C+T}}\f$
     *  \param constants Constants of the fluid.
     *
     */
    explicit VogelsLiquidDynamicViscosity(const VogelsConstants& constants)
        : _constants(constants)
    {
    }

    /// Get model name.
    std::string getName() const override
    {
        return "Liquid viscosity by Vogel's equation";
    }

    /** Get viscosity value (in SI unit: Pa * s).
     *  \param var_vals Variable values in an array. The order of its elements
     *                  is given in enum class PropertyVariableType.
     */
    double getValue(const ArrayType& var_vals) const override
    {
        const double T = var_vals[static_cast<int>(PropertyVariableType::T)];
        // Note: the constant of 1.e-3 is for the SI unit conversion.
        return 1.e-3 *
               std::exp(_constants.A + _constants.B / (_constants.C + T));
    }

    /** Get the partial differential of the viscosity with respect to
     *  temperature.
     *  \param var_vals  Variable values  in an array. The order of its elements
     *                   is given in enum class PropertyVariableType.
     *  \param var       Variable type.
     */
    double getdValue(const ArrayType& var_vals,
                     const PropertyVariableType var) const override
    {
        (void)var;
        const double T = var_vals[static_cast<int>(PropertyVariableType::T)];
        const double f_buff = _constants.B / (_constants.C + T);
        // Note: the constant of 1.e-3 is for the SI unit conversion.
        return -1.e-3 * f_buff * std::exp(_constants.A + f_buff) /
               (_constants.C + T);
    }

private:
    const VogelsConstants _constants;
};

/**  Parameters A, B, C.
 *  <a href="ddbst"> http://ddbonline.ddbst.de/VogelCalculation/VogelCalculationCGI.exe</a>
 * */
struct VogelsViscosityConstantsWater
{
    VogelsViscosityConstantsWater() {}
    const double A = -3.7188;
    const double B = 578.919;
    const double C = -137.546;
};

struct VogelsViscosityConstantsCO2
{
    VogelsViscosityConstantsCO2() {}
    const double A = -24.0592;
    const double B = 28535.2;
    const double C = 1037.41;
};

struct VogelsViscosityConstantsCH4
{
    VogelsViscosityConstantsCH4() {}
    const double A = -25.5947;
    const double B = 25392;
    const double C = 969.306;
};

}  // end namespace
}  // end namespace
#endif
