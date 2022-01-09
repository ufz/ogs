/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include "GroutParameters.h"
#include "Pipe.h"
#include "RefrigerantProperties.h"

namespace ProcessLib
{
namespace HeatTransportBHE
{
namespace BHE
{
struct AdvectiveThermalResistanceCoaxial
{
    double const inner_pipe_coaxial;
    double const a_annulus;
    double const b_annulus;
};

struct PipeWallThermalResistanceCoaxial
{
    double const inner_pipe_coaxial;
    double const annulus;
};

struct GroutAndGroutSoilExchangeThermalResistanceCoaxial
{
    double const conductive_b;
    double const grout_soil;
};

inline AdvectiveThermalResistanceCoaxial calculateAdvectiveThermalResistance(
    Pipe const& inner_pipe, Pipe const& outer_pipe,
    RefrigerantProperties const& fluid, double const Nu_inner_pipe,
    double const Nu_annulus)
{
    double const hydraulic_diameter =
        coaxialPipesAnnulusDiameter(inner_pipe, outer_pipe);

    auto advective_thermal_resistance = [&](double Nu, double diameter_ratio) {
        constexpr double pi = boost::math::constants::pi<double>();
        return 1.0 / (Nu * fluid.thermal_conductivity * pi) * diameter_ratio;
    };
    return {advective_thermal_resistance(Nu_inner_pipe, 1.),
            advective_thermal_resistance(
                Nu_annulus, hydraulic_diameter / inner_pipe.outsideDiameter()),
            advective_thermal_resistance(
                Nu_annulus, hydraulic_diameter / outer_pipe.diameter)};
}

inline PipeWallThermalResistanceCoaxial calculatePipeWallThermalResistance(
    Pipe const& inner_pipe, Pipe const& outer_pipe)
{
    return {inner_pipe.wallThermalResistance(),
            outer_pipe.wallThermalResistance()};
}

inline GroutAndGroutSoilExchangeThermalResistanceCoaxial
calculateGroutAndGroutSoilExchangeThermalResistance(
    Pipe const& outer_pipe, GroutParameters const& grout_parameters,
    double const borehole_diameter)
{
    constexpr double pi = boost::math::constants::pi<double>();

    double const outer_pipe_outside_diameter = outer_pipe.outsideDiameter();
    double const chi =
        std::log(std::sqrt(borehole_diameter * borehole_diameter +
                           outer_pipe_outside_diameter *
                               outer_pipe_outside_diameter) /
                 std::sqrt(2) / outer_pipe_outside_diameter) /
        std::log(borehole_diameter / outer_pipe_outside_diameter);
    double const R_g =
        std::log(borehole_diameter / outer_pipe_outside_diameter) / 2 /
        (pi * grout_parameters.lambda_g);
    double const conductive_b = chi * R_g;
    double const grout_soil = (1 - chi) * R_g;
    return {conductive_b, grout_soil};
}
}  // namespace BHE
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
