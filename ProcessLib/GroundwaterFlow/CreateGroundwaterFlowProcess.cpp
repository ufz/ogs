/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateGroundwaterFlowProcess.h"

#include "BaseLib/FileTools.h"
#include "ProcessLib/Utils/ParseSecondaryVariables.h"
#include "ProcessLib/Utils/ProcessUtils.h"
#include "ProcessLib/CalculateSurfaceFlux/ParseCalculateSurfaceFluxData.h"
#include "GroundwaterFlowProcess.h"
#include "GroundwaterFlowProcessData.h"

#include "MeshLib/IO/readMeshFromFile.h"

namespace ProcessLib
{
namespace GroundwaterFlow
{
std::unique_ptr<Process> createGroundwaterFlowProcess(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config,
    std::string const& project_directory,
    std::string const& output_directory)
{
    //! \ogs_file_param{prj__processes__process__type}
    config.checkConfigParameter("type", "GROUNDWATER_FLOW");

    DBUG("Create GroundwaterFlowProcess.");

    // Process variable.

    //! \ogs_file_param{prj__processes__process__GROUNDWATER_FLOW__process_variables}
    auto const pv_config = config.getConfigSubtree("process_variables");

    auto process_variables = findProcessVariables(
        variables, pv_config,
        {//! \ogs_file_param_special{prj__processes__process__GROUNDWATER_FLOW__process_variables__process_variable}
         "process_variable"});

    // Hydraulic conductivity parameter.
    auto& hydraulic_conductivity = findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__GROUNDWATER_FLOW__hydraulic_conductivity}
        "hydraulic_conductivity",
        parameters, 1);

    DBUG("Use \'%s\' as hydraulic conductivity parameter.",
         hydraulic_conductivity.name.c_str());

    GroundwaterFlowProcessData process_data{hydraulic_conductivity};

    SecondaryVariableCollection secondary_variables;

    NumLib::NamedFunctionCaller named_function_caller(
        {"GWFlow_pressure"});

    ProcessLib::parseSecondaryVariables(config, secondary_variables,
                                        named_function_caller);

    std::string mesh_name;
    std::string balance_pv_name;
    std::string balance_out_fname;
    std::unique_ptr<MeshLib::Mesh> surface_mesh;
    ProcessLib::parseCalculateSurfaceFluxData(config, mesh_name, balance_pv_name,
                                        balance_out_fname);

    if (! mesh_name.empty()) // balance is optional
    {
        mesh_name = BaseLib::copyPathToFileName(mesh_name, project_directory);

        balance_out_fname =
            BaseLib::copyPathToFileName(balance_out_fname, output_directory);

        surface_mesh.reset(MeshLib::IO::readMeshFromFile(mesh_name));

        DBUG(
            "read balance meta data:\n\tbalance mesh:\"%s\"\n\tproperty name: "
            "\"%s\"\n\toutput to: \"%s\"",
            mesh_name.c_str(), balance_pv_name.c_str(),
            balance_out_fname.c_str());

        // Surface mesh and bulk mesh must have equal axial symmetry flags!
        surface_mesh->setAxiallySymmetric(mesh.isAxiallySymmetric());
    }

    return std::make_unique<GroundwaterFlowProcess>(
        mesh, std::move(jacobian_assembler), parameters, integration_order,
        std::move(process_variables), std::move(process_data),
        std::move(secondary_variables), std::move(named_function_caller),
        surface_mesh.release(), std::move(balance_pv_name),
        std::move(balance_out_fname));
}

}  // namespace GroundwaterFlow
}  // namespace ProcessLib
