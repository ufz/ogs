// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <numbers>

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
struct BoreholeGeometry
{
    /**
     * length/depth of the BHE
     * unit is m
     */
    double const length;

    /**
     * diameter of the BHE
     * unit is m
     */
    double const diameter;

    double area() const { return std::numbers::pi * diameter * diameter / 4; }
};

BoreholeGeometry createBoreholeGeometry(BaseLib::ConfigTree const& config);

}  // namespace BHE
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
