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
    auto const borehole_geometry =
        createBoreholeGeometry(config.getConfigSubtree("borehole"));

    //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__pipes}
    auto const& pipes_config = config.getConfigSubtree("pipes");
    //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__pipes__outer}
    Pipe const outer_pipe = createPipe(pipes_config.getConfigSubtree("outer"));
    Pipe const inner_pipe =
        //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__pipes__inner}
        createPipe(pipes_config.getConfigSubtree("inner"));
    const double pipe_longitudinal_dispersion_length =
        //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__pipes__longitudinal_dispersion_length}
        pipes_config.getConfigParameter<double>(
            "longitudinal_dispersion_length");
    PipeConfigurationCoaxial const pipes{inner_pipe, outer_pipe,
                                         pipe_longitudinal_dispersion_length};

    //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__grout}
    auto const grout = createGroutParameters(config.getConfigSubtree("grout"));

    //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__refrigerant}
    auto const refrigerant =
        createRefrigerantProperties(config.getConfigSubtree("refrigerant"));

    auto const flowAndTemperatureControl = createFlowAndTemperatureControl(
        //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__flow_and_temperature_control}
        config.getConfigSubtree("flow_and_temperature_control"),
        curves,
        refrigerant);

    return {borehole_geometry, refrigerant, grout, flowAndTemperatureControl,
            pipes};
}
}  // namespace BHE
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
