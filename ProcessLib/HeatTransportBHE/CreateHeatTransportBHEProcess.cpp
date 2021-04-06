/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateHeatTransportBHEProcess.h"

#include <vector>

#include "BHE/BHETypes.h"
#include "BHE/CreateBHE1PType.h"
#include "BHE/CreateBHECoaxial.h"
#include "BHE/CreateBHEUType.h"
#include "HeatTransportBHEProcess.h"
#include "HeatTransportBHEProcessData.h"
#include "ParameterLib/Utils.h"
#include "ProcessLib/Output/CreateSecondaryVariables.h"
#ifdef OGS_USE_PYTHON
#include <pybind11/pybind11.h>
#endif  // OGS_USE_PYTHON

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
    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const& media)
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
                "Found a process variable name '{:s}'. It should be "
                "'temperature_soil' or 'temperature_BHE_X'");
        }
        auto variable = std::find_if(variables.cbegin(), variables.cend(),
                                     [&pv_name](ProcessVariable const& v) {
                                         return v.getName() == pv_name;
                                     });

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

        if (bhe_type == "1P")
        {
            bhes.emplace_back(
                BHE::createBHE1PType<BHE::BHE_1P>(bhe_config, curves));
            continue;
        }
        OGS_FATAL("Unknown BHE type '{:s}'.", bhe_type);
    }
    // end of reading BHE parameters -------------------------------------------

    auto media_map =
        MaterialPropertyLib::createMaterialSpatialDistributionMap(media, mesh);

    // find if bhe uses python boundary condition
    auto const using_python_bcs =
        visit([](auto const& bhe) { return bhe.use_python_bcs; }, bhes[0]);

    //! Python object computing BC values.
    BHEInflowPythonBoundaryConditionPythonSideInterface* py_object = nullptr;
    // create a pythonBoundaryCondition object
    if (using_python_bcs)
    {
#ifdef OGS_USE_PYTHON
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
#else
        OGS_FATAL(
            "Input files suggest the coupling of BHE with pipe network. "
            "This means the compiling flag OGS_USE_PYTHON must be switched "
            "on. ");
#endif  // OGS_USE_PYTHON
    }

    HeatTransportBHEProcessData process_data(
        std::move(media_map), std::move(bhes), py_object);

    SecondaryVariableCollection secondary_variables;

    ProcessLib::createSecondaryVariables(config, secondary_variables);

    return std::make_unique<HeatTransportBHEProcess>(
        std::move(name), mesh, std::move(jacobian_assembler), parameters,
        integration_order, std::move(process_variables),
        std::move(process_data), std::move(secondary_variables));
}
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
