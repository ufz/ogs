/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateTESProcess.h"
#include "ProcessLib/Output/CreateSecondaryVariables.h"
#include "ProcessLib/Utils/ProcessUtils.h"
#include "TESProcess.h"

namespace ProcessLib
{
namespace TES
{
std::unique_ptr<Process> createTESProcess(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{prj__processes__process__type}
    config.checkConfigParameter("type", "TES");

    DBUG("Create TESProcess.");

    //! \ogs_file_param{prj__processes__process__TES__process_variables}
    auto const pv_config = config.getConfigSubtree("process_variables");

    auto per_process_variables = findProcessVariables(
        variables, pv_config,
        {
        //! \ogs_file_param_special{prj__processes__process__TES__process_variables__fluid_pressure}
        "fluid_pressure",
        //! \ogs_file_param_special{prj__processes__process__TES__process_variables__temperature}
        "temperature",
        //! \ogs_file_param_special{prj__processes__process__TES__process_variables__vapour_mass_fraction}
        "vapour_mass_fraction"});
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>
        process_variables;
    process_variables.push_back(std::move(per_process_variables));

    SecondaryVariableCollection secondary_variables;

    NumLib::NamedFunctionCaller named_function_caller(
        {"TES_pressure", "TES_temperature", "TES_vapour_mass_fraction"});

    ProcessLib::createSecondaryVariables(config, secondary_variables,
                                         named_function_caller);

    return std::make_unique<TESProcess>(
        mesh, std::move(jacobian_assembler), parameters, integration_order,
        std::move(process_variables), std::move(secondary_variables),
        std::move(named_function_caller), config);
}

}  // namespace TES
}  // namespace ProcessLib
