// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

namespace BaseLib
{
class ConfigTree;
}
namespace ProcessLib
{
namespace HeatTransportBHE
{
namespace BHE
{
struct RefrigerantProperties
{
    /**
     * unit is kg m-1 sec-1
     */
    double const dynamic_viscosity;

    /**
     * unit is kg m-3
     */
    double const density;

    /**
     * unit is kg m sec^-3 K^-1
     */
    double const thermal_conductivity;

    /**
     * unit is m^2 sec^-2 K^-1
     */
    double const specific_heat_capacity;

    double const reference_temperature;
};

RefrigerantProperties createRefrigerantProperties(
    BaseLib::ConfigTree const& config);

}  // namespace BHE
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
