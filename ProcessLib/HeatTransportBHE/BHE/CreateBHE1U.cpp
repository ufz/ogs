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

#include "CreateBHE1U.h"
#include "BaseLib/ConfigTree.h"
#include "CreateFlowAndTemperatureControl.h"
namespace ProcessLib
{
namespace HeatTransportBHE
{
namespace BHE
{
BHE::BHE_1U createBHE1U(
    BaseLib::ConfigTree const& config,
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
        curves)
{
    auto const& borehole_config = config.getConfigSubtree("borehole");
    const double borehole_length =
        borehole_config.getConfigParameter<double>("length");
    const double borehole_diameter =
        borehole_config.getConfigParameter<double>("diameter");
    BoreholeGeometry const borehole{borehole_length, borehole_diameter};

    auto const& pipes_config = config.getConfigSubtree("pipes");
    Pipe const inlet_pipe = createPipe(pipes_config.getConfigSubtree("inlet"));
    Pipe const outlet_pipe =
        createPipe(pipes_config.getConfigSubtree("outlet"));
    const double pipe_distance =
        pipes_config.getConfigParameter<double>("distance_between_pipes");
    const double pipe_longitudinal_dispersion_length =
        pipes_config.getConfigParameter<double>(
            "longitudinal_dispersion_length");
    PipeConfiguration1U const pipes{inlet_pipe, outlet_pipe, pipe_distance,
                                    pipe_longitudinal_dispersion_length};

    auto const& grout_config = config.getConfigSubtree("grout");
    const double grout_density =
        grout_config.getConfigParameter<double>("density");
    const double grout_porosity =
        grout_config.getConfigParameter<double>("porosity");
    const double grout_heat_capacity =
        grout_config.getConfigParameter<double>("heat_capacity");
    const double grout_thermal_conductivity =
        grout_config.getConfigParameter<double>("thermal_conductivity");
    GroutParameters const grout{grout_density, grout_porosity,
                                grout_heat_capacity,
                                grout_thermal_conductivity};

    auto const& refrigerant_config = config.getConfigSubtree("refrigerant");
    double const refrigerant_density =
        refrigerant_config.getConfigParameter<double>("density");
    double const refrigerant_viscosity =
        refrigerant_config.getConfigParameter<double>("viscosity");
    double const refrigerant_heat_capacity =
        refrigerant_config.getConfigParameter<double>("specific_heat_capacity");
    double const refrigerant_thermal_conductivity =
        refrigerant_config.getConfigParameter<double>("thermal_conductivity");
    double const refrigerant_reference_temperature =
        refrigerant_config.getConfigParameter<double>("reference_temperature");
    RefrigerantProperties const refrigerant{
        refrigerant_viscosity, refrigerant_density,
        refrigerant_thermal_conductivity, refrigerant_heat_capacity,
        refrigerant_reference_temperature};

    auto const flowAndTemperatureControl = createFlowAndTemperatureControl(
        config.getConfigSubtree("flow_and_temperature_control"),
        curves,
        refrigerant);

    return {borehole, refrigerant, grout, flowAndTemperatureControl, pipes};
}
}  // namespace BHE
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
