/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

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
    double const diameter;
    double const wall_thickness;
    double const wall_thermal_conductivity;

    /// Area of the pipe's inside without the wall.
    double area() const { return circleArea(diameter); }

    /// Area of the pipe's outside including the wall thickness.
    double outsideArea() const { return circleArea(outsideDiameter()); }

    double outsideDiameter() const { return diameter + 2 * wall_thickness; }

    double wallThermalResistance() const
    {
        double const outside_diameter = outsideDiameter();

        return std::log(outside_diameter / diameter) /
               (2.0 * std::numbers::pi * wall_thermal_conductivity);
    }

private:
    double circleArea(double const diameter) const
    {
        return std::numbers::pi * diameter * diameter / 4;
    }
};

inline double coaxialPipesAnnulusDiameter(Pipe const& inner_pipe,
                                          Pipe const& outer_pipe)
{
    return outer_pipe.diameter - inner_pipe.diameter -
           2 * inner_pipe.wall_thickness;
}

Pipe createPipe(BaseLib::ConfigTree const& config);
}  // namespace BHE
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
