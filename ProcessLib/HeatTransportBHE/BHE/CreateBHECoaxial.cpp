/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "CreateBHECoaxial.h"

#include "BHE_CXA.h"
#include "BHE_CXC.h"
#include "BaseLib/ConfigTree.h"
#include "CreateFlowAndTemperatureControl.h"
namespace ProcessLib
{
namespace HeatTransportBHE
{
namespace BHE
{
static std::tuple<BoreholeGeometry,
                  RefrigerantProperties,
                  GroutParameters,
                  FlowAndTemperatureControl,
                  PipeConfigurationCoaxial,
                  bool>
parseBHECoaxialConfig(
    BaseLib::ConfigTree const& config,
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
        curves)
{
    // if the BHE is using python boundary condition
    auto const bhe_if_use_python_bc_conf =
        //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__use_bhe_pipe_network}
        config.getConfigParameter<bool>("use_bhe_pipe_network", false);
    DBUG("If using python boundary condition : {:s}",
         (bhe_if_use_python_bc_conf) ? "true" : "false");

    auto const borehole_geometry =
        //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__borehole}
        createBoreholeGeometry(config.getConfigSubtree("borehole"));

    //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__pipes}
    auto const& pipes_config = config.getConfigSubtree("pipes");
    //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__pipes__outer}
    Pipe const outer_pipe = createPipe(pipes_config.getConfigSubtree("outer"));
    Pipe const inner_pipe =
        //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__pipes__inner}
        createPipe(pipes_config.getConfigSubtree("inner"));
    const auto pipe_longitudinal_dispersion_length =
        //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__pipes__longitudinal_dispersion_length}
        pipes_config.getConfigParameter<double>(
            "longitudinal_dispersion_length");
    PipeConfigurationCoaxial const pipes{inner_pipe, outer_pipe,
                                         pipe_longitudinal_dispersion_length};

    //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__grout}
    auto const grout = createGroutParameters(config.getConfigSubtree("grout"));

    auto const refrigerant =
        //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__refrigerant}
        createRefrigerantProperties(config.getConfigSubtree("refrigerant"));

    auto const flowAndTemperatureControl = createFlowAndTemperatureControl(
        //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__flow_and_temperature_control}
        config.getConfigSubtree("flow_and_temperature_control"),
        curves,
        refrigerant);

    return {borehole_geometry,         refrigerant, grout,
            flowAndTemperatureControl, pipes,       bhe_if_use_python_bc_conf};
}

template <typename T_BHE>
T_BHE createBHECoaxial(
    BaseLib::ConfigTree const& config,
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
        curves)
{
    auto coaxial = parseBHECoaxialConfig(config, curves);
    return {std::get<0>(coaxial), std::get<1>(coaxial), std::get<2>(coaxial),
            std::get<3>(coaxial), std::get<4>(coaxial), std::get<5>(coaxial)};
}

template BHE_CXA createBHECoaxial<BHE_CXA>(
    BaseLib::ConfigTree const& config,
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
        curves);

template BHE_CXC createBHECoaxial<BHE_CXC>(
    BaseLib::ConfigTree const& config,
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
        curves);
}  // namespace BHE
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
