/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file   createLiquidFlowProcess.cpp
 *
 * Created on August 19, 2016, 1:30 PM
 */
#include "createLiquidFlowProcess.h"

#include "ProcessLib/Utils/ParseSecondaryVariables.h"
#include "ProcessLib/Utils/ProcessUtils.h"

#include "LiquidFlowProcess.h"

namespace ProcessLib
{
namespace LiquidFlow
{
std::unique_ptr<Process> createLiquidFlowProcess(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{process__type}
    config.checkConfigParameter("type", "LIQUID_FLOW");

    DBUG("Create LiquidFlowProcess.");

    // Process variable.
    auto process_variables = findProcessVariables(
        variables, config,
        {//! \ogs_file_param_special{process__LIQUID_FLOW__process_variables__process_variable}
         "process_variable"});

    SecondaryVariableCollection secondary_variables;

    NumLib::NamedFunctionCaller named_function_caller({"LiquidFlow_pressure"});

    ProcessLib::parseSecondaryVariables(config, secondary_variables,
                                        named_function_caller);

    auto const gravitational_term =
        //! \ogs_file_param{process__LIQUID_FLOW__gravitational_term}
        config.getConfigParameter<std::string>("gravitational_term");
    const bool has_gravitational_term = (gravitational_term == "enabled");

    //! \ogs_file_param{process__LIQUID_FLOW__material_property}
    auto const& mat_config = config.getConfigSubtree("material_property");

    return std::unique_ptr<Process>{new LiquidFlowProcess{
        mesh, std::move(jacobian_assembler), parameters,
        integration_order,
        std::move(process_variables), std::move(secondary_variables),
        std::move(named_function_caller), has_gravitational_term, mat_config}};
}

}  // end of namespace
}  // end of namespace
