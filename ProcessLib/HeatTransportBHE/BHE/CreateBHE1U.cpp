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
    //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger}
    auto const& borehole_heat_exchange_config = config.getConfigSubtree("borehole_heat_exchange");
    // if the BHE is using python boundary condition
    bool bhe_if_use_python_bc = false;
    if (auto const bhe_if_use_python_bc_conf =
            borehole_heat_exchange_config.getConfigParameterOptional<bool>(
                "bhe_if_use_python_bc"))
    {
        DBUG("If  using python boundary condition : %s",
             (*bhe_if_use_python_bc_conf) ? "true" : "false");
        bhe_if_use_python_bc = *bhe_if_use_python_bc_conf;
    }

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
    //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__pipes__inlet}
    Pipe const inlet_pipe = createPipe(pipes_config.getConfigSubtree("inlet"));
    Pipe const outlet_pipe =
        //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__pipes__outlet}
        createPipe(pipes_config.getConfigSubtree("outlet"));
    const double pipe_distance =
        //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__pipes__distance_between_pipes}
        pipes_config.getConfigParameter<double>("distance_between_pipes");
    const double pipe_longitudinal_dispersion_length =
        //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__pipes__longitudinal_dispersion_length}
        pipes_config.getConfigParameter<double>(
            "longitudinal_dispersion_length");
    PipeConfiguration1U const pipes{inlet_pipe, outlet_pipe, pipe_distance,
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

    return {borehole, refrigerant, grout, flowAndTemperatureControl, pipes, bhe_if_use_python_bc};
}
}  // namespace BHE
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
