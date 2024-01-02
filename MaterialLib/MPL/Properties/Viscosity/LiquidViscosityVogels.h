/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on March 4, 2021, 3:05 PM
 */

#pragma once

#include "MaterialLib/MPL/Property.h"

namespace MaterialPropertyLib
{
class Phase;

/** A temperature dependent viscosity model.
 * <a
 * href="http://ddbonline.ddbst.de/VogelCalculation/VogelCalculationCGI.exe">ddbst</a>
 */
template <typename VogelsConstants>
class LiquidViscosityVogels final : public Property
{
public:
    explicit LiquidViscosityVogels(std::string name, VogelsConstants constants)
        : constants_(std::move(constants))
    {
        name_ = std::move(name);
    };

    void checkScale() const override
    {
        if (!std::holds_alternative<Phase*>(scale_))
        {
            OGS_FATAL(
                "The property 'LiquidViscosityVogels' is "
                "implemented on the 'Phase' scale only.");
        }
    }

    /// \return The liquid viscosity.
    PropertyDataType value(VariableArray const& variable_array,
                           ParameterLib::SpatialPosition const& pos,
                           double const t, double const dt) const override;
    /// \return The derivative of the liquid viscosity with respect to
    /// temperature or phase (water) pressure.
    PropertyDataType dValue(VariableArray const& variable_array,
                            Variable const variable,
                            ParameterLib::SpatialPosition const& pos,
                            double const t, double const dt) const override;

private:
    const VogelsConstants constants_;
    // Coefficients Hi and Hij are given in two static arrays in the cpp file.
};

/**  Parameters A, B, C.
 *  <a
 * href="http://ddbonline.ddbst.de/VogelCalculation/VogelCalculationCGI.exe">ddbst</a>
 * */
struct VogelsViscosityConstantsWater
{
    VogelsViscosityConstantsWater() = default;
    const double A = -3.7188;
    const double B = 578.919;
    const double C = -137.546;
};

struct VogelsViscosityConstantsCO2
{
    VogelsViscosityConstantsCO2() = default;
    const double A = -24.0592;
    const double B = 28535.2;
    const double C = 1037.41;
};

struct VogelsViscosityConstantsCH4
{
    VogelsViscosityConstantsCH4() = default;
    const double A = -25.5947;
    const double B = 25392;
    const double C = 969.306;
};

extern template class LiquidViscosityVogels<VogelsViscosityConstantsWater>;
extern template class LiquidViscosityVogels<VogelsViscosityConstantsCO2>;
extern template class LiquidViscosityVogels<VogelsViscosityConstantsCH4>;
}  // namespace MaterialPropertyLib
