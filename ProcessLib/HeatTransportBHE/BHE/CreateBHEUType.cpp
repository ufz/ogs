// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "CreateBHEUType.h"

#include <algorithm>
#include <cmath>

#include "BHE_1U.h"
#include "BHE_2U.h"
#include "BaseLib/ConfigTree.h"
#include "BaseLib/Error.h"
#include "CreateFlowAndTemperatureControl.h"
namespace ProcessLib
{
namespace HeatTransportBHE
{
namespace BHE
{
static std::tuple<BoreholeGeometry, RefrigerantProperties, GroutParameters,
                  FlowAndTemperatureControl, PipeConfigurationUType, bool>
parseBHEUTypeConfig(
    BaseLib::ConfigTree const& config,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>>& parameters,
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
        curves,
    std::vector<MeshLib::Node*> const& bhe_nodes)
{
    // if the BHE is using python boundary condition
    auto const bhe_if_use_python_bc_conf =
        //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__use_bhe_pipe_network}
        config.getConfigParameter<bool>("use_bhe_pipe_network", false);
    DBUG("If using python boundary condition : {:s}",
         (bhe_if_use_python_bc_conf) ? "true" : "false");

    auto const borehole_geometry =
        //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__borehole}
        createBoreholeGeometry(config.getConfigSubtree("borehole"), parameters,
                               bhe_nodes);

    //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__pipes}
    auto const& pipes_config = config.getConfigSubtree("pipes");
    //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__pipes__inlet}
    Pipe const inlet_pipe = createPipe(pipes_config.getConfigSubtree("inlet"));
    Pipe const outlet_pipe =
        //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__pipes__outlet}
        createPipe(pipes_config.getConfigSubtree("outlet"));
    const auto pipe_distance =
        //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__pipes__distance_between_pipes}
        pipes_config.getConfigParameter<double>("distance_between_pipes");
    const auto pipe_longitudinal_dispersion_length =
        //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__pipes__longitudinal_dispersion_length}
        pipes_config.getConfigParameter<double>(
            "longitudinal_dispersion_length");
    PipeConfigurationUType const pipes{inlet_pipe, outlet_pipe, pipe_distance,
                                       pipe_longitudinal_dispersion_length};

    if (pipe_distance <= 0)
    {
        OGS_FATAL(
            "distance_between_pipes must be positive for U-type BHEs, got "
            "{:g}.",
            pipe_distance);
    }

    double const d0 =
        std::max(pipes.inlet.outsideDiameter(), pipes.outlet.outsideDiameter());
    if (pipe_distance < d0)
    {
        OGS_FATAL(
            "distance_between_pipes ({:g}) must be >= pipe outside diameter "
            "({:g}) for valid U-type thermal resistance formulas.",
            pipe_distance, d0);
    }

    for (int section_index = 0;
         section_index < borehole_geometry.sections.getNumberOfSections();
         ++section_index)
    {
        double const D =
            borehole_geometry.sections.diameterAtSection(section_index);
        if (D <= 2.0 * d0)
        {
            OGS_FATAL(
                "Invalid U-type geometry at section {:d}: borehole diameter "
                "{:g} must be greater than 2*pipe outside diameter {:g}.",
                section_index, D, 2.0 * d0);
        }

        double const acosh_argument =
            (D * D + d0 * d0 - pipe_distance * pipe_distance) / (2.0 * D * d0);
        if (!std::isfinite(acosh_argument) || acosh_argument < 1.0)
        {
            OGS_FATAL(
                "Invalid U-type geometry at section {:d}: acosh argument "
                "for grout resistance is {:g}, must be >= 1. "
                "(D={:g}, d0={:g}, distance_between_pipes={:g})",
                section_index, acosh_argument, D, d0, pipe_distance);
        }
    }

    //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__grout}
    auto const grout = createGroutParameters(config.getConfigSubtree("grout"));

    auto const refrigerant =
        //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__refrigerant}
        createRefrigerantProperties(config.getConfigSubtree("refrigerant"));

    auto const flowAndTemperatureControl = createFlowAndTemperatureControl(
        //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__flow_and_temperature_control}
        config.getConfigSubtree("flow_and_temperature_control"),
        parameters,
        curves,
        refrigerant);

    return {borehole_geometry,         refrigerant, grout,
            flowAndTemperatureControl, pipes,       bhe_if_use_python_bc_conf};
}

template <typename T_BHE>
T_BHE createBHEUType(
    BaseLib::ConfigTree const& config,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>>& parameters,
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
        curves,
    std::vector<MeshLib::Node*> const& bhe_nodes)
{
    auto UType = parseBHEUTypeConfig(config, parameters, curves, bhe_nodes);
    return {std::get<0>(UType), std::get<1>(UType), std::get<2>(UType),
            std::get<3>(UType), std::get<4>(UType), std::get<5>(UType)};
}

template BHE_1U createBHEUType<BHE_1U>(
    BaseLib::ConfigTree const& config,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>>& parameters,
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
        curves,
    std::vector<MeshLib::Node*> const& bhe_nodes);

template BHE_2U createBHEUType<BHE_2U>(
    BaseLib::ConfigTree const& config,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>>& parameters,
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
        curves,
    std::vector<MeshLib::Node*> const& bhe_nodes);
}  // namespace BHE
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
