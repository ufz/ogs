/**
 * \file
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateHeatTransportBHEProcess.h"

#include <vector>

#include "ParameterLib/Utils.h"
#include "ProcessLib/Output/CreateSecondaryVariables.h"

#include "BHE/BHETypes.h"
#include "BHE/CreateBHECoaxial.h"
#include "BHE/CreateBHEUType.h"
#include "HeatTransportBHEProcess.h"
#include "HeatTransportBHEProcessData.h"

namespace ProcessLib
{
namespace HeatTransportBHE
{
std::unique_ptr<Process> createHeatTransportBHEProcess(
    std::string name,
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config,
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
        curves,
    std::map<int, std::unique_ptr<MaterialPropertyLib::Medium>> const& media)
{
    //! \ogs_file_param{prj__processes__process__type}
    config.checkConfigParameter("type", "HEAT_TRANSPORT_BHE");

    DBUG("Create HeatTransportBHE Process.");

    // Process variable.

    //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__process_variables}
    auto const pv_config = config.getConfigSubtree("process_variables");
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>
        process_variables;

    // reading primary variables for each
    // BHE----------------------------------------------------------
    auto range =
        //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__process_variables__process_variable}
        pv_config.getConfigParameterList<std::string>("process_variable");
    std::vector<std::reference_wrapper<ProcessVariable>> per_process_variables;

    for (std::string const& pv_name : range)
    {
        if (pv_name != "temperature_soil" &&
            pv_name.find("temperature_BHE") == std::string::npos)
        {
            OGS_FATAL(
                "Found a process variable name '%s'. It should be "
                "'temperature_soil' or 'temperature_BHE_X'");
        }
        auto variable = std::find_if(variables.cbegin(), variables.cend(),
                                     [&pv_name](ProcessVariable const& v) {
                                         return v.getName() == pv_name;
                                     });

        if (variable == variables.end())
        {
            OGS_FATAL(
                "Could not find process variable '%s' in the provided "
                "variables "
                "list for config tag <%s>.",
                pv_name.c_str(), "process_variable");
        }
        DBUG("Found process variable '%s' for config tag <%s>.",
             variable->getName().c_str(), "process_variable");

        per_process_variables.emplace_back(
            const_cast<ProcessVariable&>(*variable));
    }
    process_variables.push_back(std::move(per_process_variables));
    // end of reading primary variables for each
    // BHE----------------------------------------------------------

    // reading BHE parameters --------------------------------------------------
    std::vector<BHE::BHETypes> bhes;
    auto const& bhe_configs =
        //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers}
        config.getConfigSubtree("borehole_heat_exchangers");

    for (
        auto const& bhe_config :
        //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger}
        bhe_configs.getConfigSubtreeList("borehole_heat_exchanger"))
    {
        // read in the parameters
        const std::string bhe_type =
            //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__type}
            bhe_config.getConfigParameter<std::string>("type");

        if (bhe_type == "1U")
        {
            bhes.emplace_back(
                BHE::createBHEUType<BHE::BHE_1U>(bhe_config, curves));
            continue;
        }

        if (bhe_type == "CXA")
        {
            bhes.emplace_back(
                BHE::createBHECoaxial<BHE::BHE_CXA>(bhe_config, curves));
            continue;
        }

        if (bhe_type == "CXC")
        {
            bhes.emplace_back(
                BHE::createBHECoaxial<BHE::BHE_CXC>(bhe_config, curves));
            continue;
        }

        if (bhe_type == "2U")
        {
            bhes.emplace_back(
                BHE::createBHEUType<BHE::BHE_2U>(bhe_config, curves));
            continue;
        }
        OGS_FATAL("Unknown BHE type '%s'.", bhe_type.c_str());
    }
    // end of reading BHE parameters -------------------------------------------

    auto media_map =
        MaterialPropertyLib::createMaterialSpatialDistributionMap(media, mesh);

    HeatTransportBHEProcessData process_data{std::move(media_map),
                                             std::move(bhes)};

    SecondaryVariableCollection secondary_variables;

    NumLib::NamedFunctionCaller named_function_caller(
        {"HeatTransportBHE_Temperature"});

    ProcessLib::createSecondaryVariables(config, secondary_variables,
                                         named_function_caller);

    return std::make_unique<HeatTransportBHEProcess>(
        std::move(name), mesh, std::move(jacobian_assembler), parameters,
        integration_order, std::move(process_variables),
        std::move(process_data), std::move(secondary_variables),
        std::move(named_function_caller));
}
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
