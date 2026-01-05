// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "CreateHeatTransportBHEProcess.h"

#include <pybind11/pybind11.h>

#include <algorithm>
#include <map>
#include <ranges>
#include <vector>

#include "BHE/BHETypes.h"
#include "BHE/CreateBHE1PType.h"
#include "BHE/CreateBHECoaxial.h"
#include "BHE/CreateBHEUType.h"
#include "HeatTransportBHEProcess.h"
#include "HeatTransportBHEProcessData.h"
#include "MaterialLib/Utils/MediaCreation.h"
#include "ParameterLib/Utils.h"
#include "ProcessLib/HeatTransportBHE/BHE/MeshUtils.h"
#include "ProcessLib/Output/CreateSecondaryVariables.h"

namespace ProcessLib
{
namespace HeatTransportBHE
{
std::map<std::string_view,
         std::function<BHE::BHETypes(
             BaseLib::ConfigTree const&,
             std::map<std::string,
                      std::unique_ptr<
                          MathLib::PiecewiseLinearInterpolation>> const&)>>
    bheCreators = {{"1U",
                    [](auto& config, auto& curves) {
                        return BHE::BHE_1U(
                            BHE::createBHEUType<BHE::BHE_1U>(config, curves));
                    }},
                   {"2U",
                    [](auto& config, auto& curves) {
                        return BHE::BHE_2U(
                            BHE::createBHEUType<BHE::BHE_2U>(config, curves));
                    }},
                   {"CXA",
                    [](auto& config, auto& curves) {
                        return BHE::BHE_CXA(BHE::createBHECoaxial<BHE::BHE_CXA>(
                            config, curves));
                    }},
                   {"CXC",
                    [](auto& config, auto& curves) {
                        return BHE::BHE_CXC(BHE::createBHECoaxial<BHE::BHE_CXC>(
                            config, curves));
                    }},
                   {"1P", [](auto& config, auto& curves) {
                        return BHE::BHE_1P(
                            BHE::createBHE1PType<BHE::BHE_1P>(config, curves));
                    }}};

void createAndInsertBHE(
    const std::string& bhe_type, const std::vector<int>& bhe_ids_of_this_bhe,
    const BaseLib::ConfigTree& bhe_config,
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
        curves,
    std::map<int, BHE::BHETypes>& bhes_map)
{
    auto bhe_creator_it = bheCreators.find(bhe_type);
    if (bhe_creator_it == bheCreators.end())
    {
        OGS_FATAL("Unknown BHE type: {:s}", bhe_type);
    }
    for (auto const& id : bhe_ids_of_this_bhe)
    {
        std::pair<std::map<int, BHE::BHETypes>::iterator, bool> result;
        if (id == bhe_ids_of_this_bhe[0])
        {
            result = bhes_map.try_emplace(
                id, bhe_creator_it->second(bhe_config, curves));
        }
        else
        {
            result = bhes_map.try_emplace(
                id, bhes_map.find(bhe_ids_of_this_bhe[0])->second);
        }
        if (!result.second)
        {
            OGS_FATAL(
                "BHE with id '{:d}' is already present in the list! Check for "
                "duplicate definitions of BHE ids.",
                id);
        }
    }
}

std::unique_ptr<Process> createHeatTransportBHEProcess(
    std::string const& name,
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config,
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
        curves,
    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const& media)
{
    //! \ogs_file_param{prj__processes__process__type}
    config.checkConfigParameter("type", "HEAT_TRANSPORT_BHE");

    DBUG("Create HeatTransportBHE Process.");

    /// \section processvariablesbhe Process Variables

    //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__process_variables}
    auto const pv_config = config.getConfigSubtree("process_variables");
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>
        process_variables;

    // reading primary variables for each
    // BHE----------------------------------------------------------
    /// Primary process variables as they appear in the global component vector:
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
                "Found a process variable name '{}'. It should be "
                "'temperature_soil' or 'temperature_BHE_X'",
                pv_name);
        }
        auto variable = std::find_if(variables.cbegin(), variables.cend(),
                                     [&pv_name](ProcessVariable const& v)
                                     { return v.getName() == pv_name; });

        if (variable == variables.end())
        {
            OGS_FATAL(
                "Could not find process variable '{:s}' in the provided "
                "variables list for config tag <{:s}>.",
                pv_name, "process_variable");
        }
        DBUG("Found process variable '{:s}' for config tag <{:s}>.",
             variable->getName(), "process_variable");

        per_process_variables.emplace_back(
            const_cast<ProcessVariable&>(*variable));
    }
    process_variables.push_back(std::move(per_process_variables));
    // end of reading primary variables for each
    // BHE----------------------------------------------------------

    /// \section parametersbhe Process Parameters
    // reading BHE parameters --------------------------------------------------

    auto bhe_mesh_data = getBHEDataInMesh(mesh);

    auto const& bhe_configs =
        //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers}
        config.getConfigSubtree("borehole_heat_exchangers");

    auto const using_server_communication =
        //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__use_server_communication}
        config.getConfigParameter<bool>("use_server_communication", false);

    auto const mass_lumping =
        //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__mass_lumping}
        config.getConfigParameter<bool>("mass_lumping", false);

    auto const using_algebraic_bc =
        //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__use_algebraic_bc}
        config.getConfigParameter<bool>("use_algebraic_bc", false);

    auto const weighting_factor =
        //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__weighting_factor}
        config.getConfigParameter<float>("weighting_factor", 1000.0);

    auto const is_linear =
        //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__linear}
        config.getConfigParameter<bool>("linear", false);
    if (is_linear)
    {
        if (!using_algebraic_bc)
        {
            WARN(
                "You specified that the process simulated by OGS is linear. "
                "With that optimization the process will be assembled only "
                "once and the non-linear solver will only do iterations per "
                "time step to fulfill the BHE boundary conditions. No other "
                "non-linearities will be resolved and OGS will not detect if "
                "there are any non-linearities. It is your responsibility to "
                "ensure that the assembled equation systems are linear, "
                "indeed! There is no safety net!");
        }
        else
        {
            WARN(
                "You specified that the process simulated by OGS is linear. "
                "With that optimization the process will be assembled only "
                "once and the non-linear solver will do only one iteration per "
                "time step. No non-linearities will be resolved and OGS will "
                "not detect if there are any non-linearities. It is your "
                "responsibility to ensure that the assembled equation systems "
                "are linear, indeed! There is no safety net!");
        }
    }

    std::map<int, BHE::BHETypes> bhes_map;

    int bhe_iterator = 0;

    for (
        auto const& bhe_config :
        //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger}
        bhe_configs.getConfigSubtreeList("borehole_heat_exchanger"))
    {
        auto bhe_id_string =
            //! \ogs_file_attr{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__id}
            bhe_config.getConfigAttribute<std::string>(
                "id", std::to_string(bhe_iterator));

        std::vector<int> bhe_ids_of_this_bhe;

        if (bhe_id_string == "*")
        {
            int size = static_cast<int>(bhe_mesh_data.BHE_mat_IDs.size());
            bhe_ids_of_this_bhe.resize(size);
            std::iota(bhe_ids_of_this_bhe.begin(), bhe_ids_of_this_bhe.end(),
                      0);
        }
        else
        {
            bhe_ids_of_this_bhe =
                MaterialLib::splitMaterialIdString(bhe_id_string);
        }

        // read in the parameters
        const std::string bhe_type =
            //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__type}
            bhe_config.getConfigParameter<std::string>("type");

        createAndInsertBHE(bhe_type, bhe_ids_of_this_bhe, bhe_config, curves,
                           bhes_map);
        bhe_iterator++;
    }

    if (static_cast<int>(bhes_map.size()) - 1 != bhes_map.rbegin()->first)
    {
        OGS_FATAL(
            "The maximum given BHE id '{:d}' did not match the number of given "
            "BHE definitions '{:d}'. The BHE ids needs to be defined starting "
            "from 0, so the maximum BHE id needs to be number of BHE "
            "definitions minus 1. After all definitions there are no gaps "
            "allowed between the given ids.",
            bhes_map.rbegin()->first, bhes_map.size());
    }

    std::vector<BHE::BHETypes> bhes;
    bhes.reserve(bhes_map.size());
    std::ranges::copy(bhes_map | std::views::values, std::back_inserter(bhes));
    //  end of reading BHE parameters
    //  -------------------------------------------

    auto media_map =
        MaterialPropertyLib::createMaterialSpatialDistributionMap(media, mesh);

    // find if bhe uses python boundary condition
    auto const using_tespy =
        visit([](auto const& bhe) { return bhe.use_python_bcs; }, bhes[0]);

    //! Python object computing BC values.
    BHEInflowPythonBoundaryConditionPythonSideInterface* py_object = nullptr;
    // create a pythonBoundaryCondition object
    if (using_tespy || using_server_communication)
    {
        // Evaluate Python code in scope of main module
        pybind11::object scope =
            pybind11::module::import("__main__").attr("__dict__");

        if (!scope.contains("bc_bhe"))
            OGS_FATAL(
                "Function 'bc_bhe' is not defined in the python script file, "
                "or there was no python script file specified.");

        py_object =
            scope["bc_bhe"]
                .cast<BHEInflowPythonBoundaryConditionPythonSideInterface*>();

        if (py_object == nullptr)
            OGS_FATAL(
                "Not able to access the correct bc pointer from python script "
                "file specified.");

        // create BHE network dataframe from Python
        py_object->dataframe_network = py_object->initializeDataContainer();
        if (!py_object->isOverriddenEssential())
        {
            DBUG(
                "Method `initializeDataContainer' not overridden in Python "
                "script.");
        }
        // clear ogs bc_node_id memory in dataframe
        std::get<3>(py_object->dataframe_network).clear();  // ogs_bc_node_id

        // here calls the tespyHydroSolver to get the pipe flow velocity in bhe
        // network
        /* for 2U type the flowrate initialization process below causes conflict
        // replace the value in flow velocity Matrix _u
        auto const tespy_flow_rate = std::get<4>(py_object->dataframe_network);
        const std::size_t n_bhe = tespy_flow_rate.size();
        if (bhes.size() != n_bhe)
            OGS_FATAL(
                "The number of BHEs defined in OGS and TESPy are not the "
                "same!");

        for (std::size_t idx_bhe = 0; idx_bhe < n_bhe; idx_bhe++)
        {
            // the flow_rate in OGS should be updated from the flow_rate
            // computed by TESPy.
            auto update_flow_rate = [&](auto& bhe) {
                bhe.updateHeatTransferCoefficients(tespy_flow_rate[idx_bhe]);
            };
            visit(update_flow_rate, bhes[idx_bhe]);
        }
        */
    }

    HeatTransportBHEProcessData process_data(
        std::move(media_map), std::move(bhes), py_object, using_tespy,
        using_server_communication, mass_lumping,
        {using_algebraic_bc, weighting_factor}, is_linear);

    SecondaryVariableCollection secondary_variables;

    ProcessLib::createSecondaryVariables(config, secondary_variables);

    return std::make_unique<HeatTransportBHEProcess>(
        std::move(name), mesh, std::move(jacobian_assembler), parameters,
        integration_order, std::move(process_variables),
        std::move(process_data), std::move(secondary_variables),
        std::move(bhe_mesh_data));
}
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
