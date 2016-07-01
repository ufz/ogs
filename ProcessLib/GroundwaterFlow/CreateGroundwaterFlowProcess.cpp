/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateGroundwaterFlowProcess.h"

#include "GroundwaterFlowProcess.h"
#include "GroundwaterFlowProcessData.h"

namespace ProcessLib
{
namespace GroundwaterFlow
{
std::unique_ptr<Process> createGroundwaterFlowProcess(
    MeshLib::Mesh& mesh,
    Process::NonlinearSolver& nonlinear_solver,
    std::unique_ptr<Process::TimeDiscretization>&& time_discretization,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{process__type}
    config.checkConfigParameter("type", "GROUNDWATER_FLOW");

    DBUG("Create GroundwaterFlowProcess.");

    // Process variable.
    auto process_variables = findProcessVariables(
        variables, config,
        {//! \ogs_file_param_special{process__GROUNDWATER_FLOW__process_variables__process_variable}
         "process_variable"});

    // Hydraulic conductivity parameter.
    auto& hydraulic_conductivity = findParameter<double,
                                                 MeshLib::Element const&>(
        config,
        //! \ogs_file_param_special{process__GROUNDWATER_FLOW__hydraulic_conductivity}
        "hydraulic_conductivity",
        parameters);

    DBUG("Use \'%s\' as hydraulic conductivity parameter.",
         hydraulic_conductivity.name.c_str());

    GroundwaterFlowProcessData process_data{hydraulic_conductivity};

    SecondaryVariableCollection secondary_variables{
        //! \ogs_file_param{process__secondary_variables}
        config.getConfigSubtreeOptional("secondary_variables"),
        {//! \ogs_file_param_special{process__GROUNDWATER_FLOW__secondary_variables__darcy_velocity_x}
         "darcy_velocity_x",
         //! \ogs_file_param_special{process__GROUNDWATER_FLOW__secondary_variables__darcy_velocity_y}
         "darcy_velocity_y",
         //! \ogs_file_param_special{process__GROUNDWATER_FLOW__secondary_variables__darcy_velocity_z}
         "darcy_velocity_z"}};

    ProcessOutput
        //! \ogs_file_param{process__output}
        process_output{config.getConfigSubtree("output"), process_variables,
                       secondary_variables};

    return std::unique_ptr<Process>{new GroundwaterFlowProcess{
        mesh, nonlinear_solver, std::move(time_discretization),
        std::move(process_variables), std::move(process_data),
        std::move(secondary_variables), std::move(process_output)}};
}

}  // namespace GroundwaterFlow
}  // namespace ProcessLib
