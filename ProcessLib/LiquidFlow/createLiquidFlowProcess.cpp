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
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{process__type}
    config.checkConfigParameter("type", "LIQUID_FLOW");

    DBUG("Create LiquidFlowProcess.");

    auto process_variables = findProcessVariables(
        //! \ogs_file_param_special{process__LIQUID_FLOW__process_variables__liquid_pressure}
        variables, config, {"liquid_pressure"});

    SecondaryVariableCollection secondary_variables;

    NumLib::NamedFunctionCaller named_function_caller({"Liquid_flow_pressure"});

    ProcessLib::parseSecondaryVariables(config, secondary_variables,
                                        named_function_caller);

    return std::unique_ptr<Process>{new LiquidFlowProcess{
        mesh, std::move(jacobian_assembler), parameters,
        std::move(process_variables), std::move(secondary_variables),
        std::move(named_function_caller), config}};
}

}  // end of namespace
}  // end of namespace
