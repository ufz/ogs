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

#include "CreateBHECXC.h"
#include "BaseLib/ConfigTree.h"
#include "CreateFlowAndTemperatureControl.h"
namespace ProcessLib
{
namespace HeatTransportBHE
{
namespace BHE
{
BHE::BHE_CXC createBHECXC(
    BaseLib::ConfigTree const& config,
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
        curves)
{
    //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__borehole}
    auto const& borehole_config = config.getConfigSubtree("borehole");
    const double borehole_length =
        //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__borehole__length}
        borehole_config.getConfigParameter<double>("length");
    const double borehole_diameter =
        //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__borehole__diameter}
        borehole_config.getConfigParameter<double>("diameter");
    BoreholeGeometry const borehole{borehole_length, borehole_diameter};

    //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__pipes}
    auto const& pipes_config = config.getConfigSubtree("pipes");
    //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__pipes__inner}
    Pipe const inner_pipe = createPipe(pipes_config.getConfigSubtree("inner"));
    Pipe const outer_pipe =
        //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__pipes__outer}
        createPipe(pipes_config.getConfigSubtree("outer"));
    const double pipe_longitudinal_dispersion_length =
        //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__pipes__longitudinal_dispersion_length}
        pipes_config.getConfigParameter<double>(
            "longitudinal_dispersion_length");
    PipeConfigurationCXC const pipes{inner_pipe, outer_pipe,
                                     pipe_longitudinal_dispersion_length};

    //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__grout}
    auto const& grout_config = config.getConfigSubtree("grout");
    const double grout_density =
        //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__grout__density}
        grout_config.getConfigParameter<double>("density");
    const double grout_porosity =
        //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__grout__porosity}
        grout_config.getConfigParameter<double>("porosity");
    const double grout_heat_capacity =
        //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__grout__heat_capacity}
        grout_config.getConfigParameter<double>("heat_capacity");
    const double grout_thermal_conductivity =
        //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__grout__thermal_conductivity}
        grout_config.getConfigParameter<double>("thermal_conductivity");
    GroutParameters const grout{grout_density, grout_porosity,
                                grout_heat_capacity,
                                grout_thermal_conductivity};

    //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__refrigerant}
    auto const& refrigerant_config = config.getConfigSubtree("refrigerant");
    double const refrigerant_density =
        //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__refrigerant__density}
        refrigerant_config.getConfigParameter<double>("density");
    double const refrigerant_viscosity =
        //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__refrigerant__viscosity}
        refrigerant_config.getConfigParameter<double>("viscosity");
    double const refrigerant_heat_capacity =
        //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__refrigerant__specific_heat_capacity}
        refrigerant_config.getConfigParameter<double>("specific_heat_capacity");
    double const refrigerant_thermal_conductivity =
        //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__refrigerant__thermal_conductivity}
        refrigerant_config.getConfigParameter<double>("thermal_conductivity");
    double const refrigerant_reference_temperature =
        //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__refrigerant__reference_temperature}
        refrigerant_config.getConfigParameter<double>("reference_temperature");
    RefrigerantProperties const refrigerant{
        refrigerant_viscosity, refrigerant_density,
        refrigerant_thermal_conductivity, refrigerant_heat_capacity,
        refrigerant_reference_temperature};

    auto const flowAndTemperatureControl = createFlowAndTemperatureControl(
        //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__flow_and_temperature_control}
        config.getConfigSubtree("flow_and_temperature_control"),
        curves,
        refrigerant);

    return {borehole, refrigerant, grout, flowAndTemperatureControl, pipes};
}
}  // namespace BHE
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
