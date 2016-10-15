/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateDensityDrivenFlowProcess.h"

#include "ProcessLib/Utils/ParseSecondaryVariables.h"
#include "ProcessLib/Utils/ProcessUtils.h"
#include "DensityDrivenFlowProcess.h"
#include "DensityDrivenFlowProcessData.h"

namespace ProcessLib
{
namespace DensityDrivenFlow
{
std::unique_ptr<Process> createDensityDrivenFlowProcess(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{process__type}
    config.checkConfigParameter("type", "DensityDrivenFlow");

    DBUG("Create DensityDrivenFlowProcess.");

    // Process variable.
    auto process_variables =
        findProcessVariables(variables, config, {
            "temperature_variable" , "pressure_variable"});  // configure two Pcs

    // Thermal conductivity parameter.
    auto& thermal_conductivity = findParameter<double>(
        config,
        //! \ogs_file_param_special{process__DensityDrivenFlow__thermal_conductivity}
        "thermal_conductivity",parameters, 1);

    DBUG("Use \'%s\' as thermal conductivity parameter.",
         thermal_conductivity.name.c_str());

    DensityDrivenFlowProcessData process_data{thermal_conductivity};

    SecondaryVariableCollection secondary_variables;

    NumLib::NamedFunctionCaller named_function_caller(
   {"DensityDrivenFlow_temperature_pressure"});

    ProcessLib::parseSecondaryVariables(config, secondary_variables,
    named_function_caller);

    return std::unique_ptr<Process>{new DensityDrivenFlowProcess{
        mesh, std::move(jacobian_assembler), parameters, integration_order,
        std::move(process_variables), std::move(process_data),
        std::move(secondary_variables), std::move(named_function_caller)}};
}

}  // namespace DensityDrivenFlow
}  // namespace ProcessLib
