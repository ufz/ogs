/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "GroutParameters.h"
#include "Physics.h"
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
    double inner_pipe_coaxial;
    double a_annulus;
    double b_annulus;
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
    static constexpr double pi = boost::math::constants::pi<double>();

    double const hydraulic_diameter_annulus = outer_pipe.diameter -
                                              inner_pipe.diameter -
                                              2 * inner_pipe.wall_thickness;
    double const inner_pipe_outside_diameter =
        PipeOutsideDiameter(inner_pipe.diameter, inner_pipe.wall_thickness);
    double const inner_pipe_coaxial =
        1.0 / (Nu_inner_pipe * fluid.thermal_conductivity * pi);
    double const a_annulus =
        1.0 / (Nu_annulus * fluid.thermal_conductivity * pi) *
        (hydraulic_diameter_annulus / inner_pipe_outside_diameter);
    double const b_annulus = 1.0 /
                             (Nu_annulus * fluid.thermal_conductivity * pi) *
                             (hydraulic_diameter_annulus / outer_pipe.diameter);
    return {inner_pipe_coaxial, a_annulus, b_annulus};
}

inline PipeWallThermalResistanceCoaxial calculatePipeWallThermalResistance(
    Pipe const& inner_pipe, Pipe const& outer_pipe)
{
    static constexpr double pi = boost::math::constants::pi<double>();

    double const inner_pipe_outside_diameter =
        PipeOutsideDiameter(inner_pipe.diameter, inner_pipe.wall_thickness);
    double const outer_pipe_outside_diameter =
        PipeOutsideDiameter(outer_pipe.diameter, outer_pipe.wall_thickness);
    double const inner_pipe_coaxial =
        std::log(inner_pipe_outside_diameter / inner_pipe.diameter) /
        (2.0 * pi * inner_pipe.wall_thermal_conductivity);
    double const Annulus =
        std::log(outer_pipe_outside_diameter / outer_pipe.diameter) /
        (2.0 * pi * outer_pipe.wall_thermal_conductivity);
    return {inner_pipe_coaxial, Annulus};
}

inline GroutAndGroutSoilExchangeThermalResistanceCoaxial
calculateGroutAndGroutSoilExchangeThermalResistance(
    Pipe const& outer_pipe, GroutParameters const& grout_parameters,
    double const borehole_diameter)
{
    static constexpr double pi = boost::math::constants::pi<double>();

    double const outer_pipe_outside_diameter =
        PipeOutsideDiameter(outer_pipe.diameter, outer_pipe.wall_thickness);
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
