// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <cmath>
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
struct Pipe
{
    /// Inner diameter [m].
    double const diameter;

    /// Wall thickness [m].
    double const wall_thickness;

    /// Wall thermal conductivity [W/(m*K)].
    double const wall_thermal_conductivity;

    double outsideDiameter() const { return diameter + 2 * wall_thickness; }

    double outsideArea() const { return circleArea(outsideDiameter()); }

    double area() const { return circleArea(diameter); }

    double wallThermalResistance() const
    {
        return std::log(outsideDiameter() / diameter) /
               (2.0 * std::numbers::pi * wall_thermal_conductivity);
    }

private:
    static double circleArea(double const d)
    {
        return std::numbers::pi * d * d / 4;
    }
};

inline double coaxialPipesAnnulusDiameter(Pipe const& inner_pipe,
                                          Pipe const& outer_pipe)
{
    return outer_pipe.diameter - inner_pipe.outsideDiameter();
}

Pipe createPipe(BaseLib::ConfigTree const& config);
}  // namespace BHE
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
