/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include "Physics.h"
#include "Pipe.h"
#include "RefrigerantProperties.h"

namespace ProcessLib
{
namespace HeatTransportBHE
{
namespace BHE
{
struct ThermoMechanicalFlowProperties
{
    double velocity;
    double nusselt_number;
};

inline ThermoMechanicalFlowProperties
calculateThermoMechanicalFlowPropertiesPipe(Pipe const& pipe,
                                            double const length,
                                            RefrigerantProperties const& fluid,
                                            double const flow_rate)
{
    double const Pr =
        prandtlNumber(fluid.dynamic_viscosity, fluid.specific_heat_capacity,
                      fluid.thermal_conductivity);

    double const velocity = flow_rate / pipe.area();
    double const Re = reynoldsNumber(velocity, pipe.diameter,
                                     fluid.dynamic_viscosity, fluid.density);
    double const nusselt_number = nusseltNumber(Re, Pr, pipe.diameter, length);
    return {velocity, nusselt_number};
}

inline ThermoMechanicalFlowProperties
calculateThermoMechanicalFlowPropertiesAnnulus(
    Pipe const& inner_pipe, Pipe const& outer_pipe, double const length,
    RefrigerantProperties const& fluid, double const flow_rate)
{
    double const Pr =
        prandtlNumber(fluid.dynamic_viscosity, fluid.specific_heat_capacity,
                      fluid.thermal_conductivity);

    double const inner_pipe_outside_diameter = inner_pipe.outsideDiameter();

    // Velocity between the outer pipe and inner pipe.
    double const velocity =
        flow_rate / (outer_pipe.area() - inner_pipe.outsideArea());

    double const Re = reynoldsNumber(
        velocity, outer_pipe.diameter - inner_pipe_outside_diameter,
        fluid.dynamic_viscosity, fluid.density);

    double const diameter_ratio =
        inner_pipe_outside_diameter / outer_pipe.diameter;
    double const pipe_aspect_ratio =
        (outer_pipe.diameter - inner_pipe_outside_diameter) / length;
    double const nusselt_number =
        nusseltNumberAnnulus(Re, Pr, diameter_ratio, pipe_aspect_ratio);
    return {velocity, nusselt_number};
}
}  // namespace BHE
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
