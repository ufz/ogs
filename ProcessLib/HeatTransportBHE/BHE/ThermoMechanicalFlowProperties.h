/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
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
    PipeParameters const& pipe, double const length,
    RefrigerantProperties const& fluid, double const flow_rate)
{
    double const Pr =
        prandtlNumber(fluid.mu_r, fluid.specific_heat_capacity, fluid.lambda_r);

    double const velocity =
        annulusFlowVelocity(flow_rate, pipe.r_outer, pipe.r_inner + pipe.b_in);

    double const Re = reynoldsNumber(
        velocity, 2.0 * (pipe.r_outer - (pipe.r_inner + pipe.b_in)), fluid.mu_r,
        fluid.density);

    double const diameter_ratio = (pipe.r_inner + pipe.b_in) / pipe.r_outer;
    double const pipe_aspect_ratio =
        (2 * pipe.r_outer - 2 * (pipe.r_inner + pipe.b_in)) / length;
    double const nusselt_number =
        nusseltNumberAnnulus(Re, Pr, diameter_ratio, pipe_aspect_ratio);
    return {velocity, nusselt_number};
}
}  // namespace BHE
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
