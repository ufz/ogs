// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "CreateBHECoaxial.h"

#include "BHE_CXA.h"
#include "BHE_CXC.h"
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
                  FlowAndTemperatureControl, PipeConfigurationCoaxial, bool>
parseBHECoaxialConfig(
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

    double const annulus_diameter =
        coaxialPipesAnnulusDiameter(pipes.inner_pipe, pipes.outer_pipe);
    if (annulus_diameter <= 0)
    {
        OGS_FATAL(
            "Invalid coaxial pipe geometry: outer pipe inner diameter ({:g}) "
            "must be greater than inner pipe outside diameter ({:g}).",
            pipes.outer_pipe.diameter, pipes.inner_pipe.outsideDiameter());
    }

    double const outer_pipe_outside_diameter =
        pipes.outer_pipe.outsideDiameter();
    for (int section_index = 0;
         section_index < borehole_geometry.sections.getNumberOfSections();
         ++section_index)
    {
        double const D =
            borehole_geometry.sections.diameterAtSection(section_index);
        if (D <= outer_pipe_outside_diameter)
        {
            OGS_FATAL(
                "Invalid coaxial geometry at section {:d}: borehole diameter "
                "{:g} must be greater than outer pipe outside diameter {:g}.",
                section_index, D, outer_pipe_outside_diameter);
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
T_BHE createBHECoaxial(
    BaseLib::ConfigTree const& config,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>>& parameters,
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
        curves,
    std::vector<MeshLib::Node*> const& bhe_nodes)
{
    auto coaxial = parseBHECoaxialConfig(config, parameters, curves, bhe_nodes);
    return {std::get<0>(coaxial), std::get<1>(coaxial), std::get<2>(coaxial),
            std::get<3>(coaxial), std::get<4>(coaxial), std::get<5>(coaxial)};
}

template BHE_CXA createBHECoaxial<BHE_CXA>(
    BaseLib::ConfigTree const& config,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>>& parameters,
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
        curves,
    std::vector<MeshLib::Node*> const& bhe_nodes);

template BHE_CXC createBHECoaxial<BHE_CXC>(
    BaseLib::ConfigTree const& config,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>>& parameters,
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
        curves,
    std::vector<MeshLib::Node*> const& bhe_nodes);
}  // namespace BHE
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
