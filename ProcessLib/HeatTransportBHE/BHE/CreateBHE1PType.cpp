/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "BHE_1P.h"
#include "BaseLib/ConfigTree.h"
#include "CreateBHEUType.h"
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
                  PipeConfiguration1PType,
                  bool>
parseBHE1PTypeConfig(
    BaseLib::ConfigTree const& config,
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
        curves)
{
    // if the BHE is using python boundary condition
    auto const bhe_if_use_python_bc_conf =
        //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__use_bhe_pipe_network}
        config.getConfigParameter<bool>("use_bhe_pipe_network", false);

    if (bhe_if_use_python_bc_conf)
    {
        DBUG("BHE 1P using python boundary conditions.");
    }

    auto const borehole_geometry =
        //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__borehole}
        createBoreholeGeometry(config.getConfigSubtree("borehole"));

    //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__pipes}
    auto const& pipes_config = config.getConfigSubtree("pipes");
    //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__pipes__inlet}
    Pipe const inlet_pipe = createPipe(pipes_config.getConfigSubtree("inlet"));

    const auto pipe_longitudinal_dispersion_length =
        //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__pipes__longitudinal_dispersion_length}
        pipes_config.getConfigParameter<double>(
            "longitudinal_dispersion_length");
    PipeConfiguration1PType const pipes{inlet_pipe,
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
T_BHE createBHE1PType(
    BaseLib::ConfigTree const& config,
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
        curves)
{
    auto SinglePipeType = parseBHE1PTypeConfig(config, curves);
    return {std::get<0>(SinglePipeType), std::get<1>(SinglePipeType),
            std::get<2>(SinglePipeType), std::get<3>(SinglePipeType),
            std::get<4>(SinglePipeType), std::get<5>(SinglePipeType)};
}

template BHE_1P createBHE1PType<BHE_1P>(
    BaseLib::ConfigTree const& config,
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
        curves);
}  // namespace BHE
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
