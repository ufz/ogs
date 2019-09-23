/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file
 *
 * Created on August 19, 2016, 1:30 PM
 */
#include "CreateLiquidFlowProcess.h"

#include <algorithm>

#include "MaterialLib/PhysicalConstant.h"
#include "MeshLib/MeshGenerators/MeshGenerator.h"
#include "ParameterLib/Utils.h"
#include "ProcessLib/Output/CreateSecondaryVariables.h"
#include "ProcessLib/Utils/ProcessUtils.h"

#include "LiquidFlowProcess.h"

namespace ProcessLib
{
namespace LiquidFlow
{
std::unique_ptr<Process> createLiquidFlowProcess(
    std::string name,
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{prj__processes__process__type}
    config.checkConfigParameter("type", "LIQUID_FLOW");

    DBUG("Create LiquidFlowProcess.");

    // Process variable.

    //! \ogs_file_param{prj__processes__process__LIQUID_FLOW__process_variables}
    auto const pv_config = config.getConfigSubtree("process_variables");

    auto per_process_variables = findProcessVariables(
        variables, pv_config,
        {//! \ogs_file_param_special{prj__processes__process__LIQUID_FLOW__process_variables__process_variable}
         "process_variable"});
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>
        process_variables;
    process_variables.push_back(std::move(per_process_variables));

    SecondaryVariableCollection secondary_variables;

    NumLib::NamedFunctionCaller named_function_caller({"LiquidFlow_pressure"});

    ProcessLib::createSecondaryVariables(config, secondary_variables,
                                         named_function_caller);

    // Get the gravity vector for the Darcy velocity
    //! \ogs_file_param{prj__processes__process__LIQUID_FLOW__darcy_gravity}
    auto const& darcy_g_config = config.getConfigSubtree("darcy_gravity");
    const auto gravity_axis_id_input =
        //! \ogs_file_param{prj__processes__process__LIQUID_FLOW__darcy_gravity__axis_id}
        darcy_g_config.getConfigParameter<int>("axis_id");
    if (gravity_axis_id_input >= static_cast<int>(mesh.getDimension()))
    {
        OGS_FATAL(
            "The gravity axis must be a number between 0 and one less than the "
            "mesh dimension, which is %d. Read gravity axis %d from input "
            "file.",
            mesh.getDimension(), gravity_axis_id_input);
    }
    const auto g =
        //! \ogs_file_param{prj__processes__process__LIQUID_FLOW__darcy_gravity__g}
        darcy_g_config.getConfigParameter<double>("g");
    if (g < 0.)
    {
        OGS_FATAL("Gravity magnitude must be non-negative.");
    }
    const int gravity_axis_id = (g == 0.) ? -1 : gravity_axis_id_input;

    auto const& refT =
        //! \ogs_file_param{prj__processes__process__LIQUID_FLOW__reference_temperature}
        config.getConfigParameterOptional<double>("reference_temperature");
    const double reference_temperature =
        refT ? *refT
             : MaterialLib::PhysicalConstant::CelsiusZeroInKelvin + 18.0;

    //! \ogs_file_param{prj__processes__process__LIQUID_FLOW__material_property}
    auto const& mat_config = config.getConfigSubtree("material_property");

    auto const material_ids = materialIDs(mesh);
    if (material_ids)
    {
        INFO("The liquid flow is in heterogeneous porous media.");
    }
    else
    {
        INFO("The liquid flow is in homogeneous porous media.");
    }

    return std::make_unique<LiquidFlowProcess>(
        std::move(name), mesh, std::move(jacobian_assembler), parameters,
        integration_order, std::move(process_variables),
        std::move(secondary_variables), std::move(named_function_caller),
        material_ids, gravity_axis_id, g, reference_temperature, mat_config);
}

}  // namespace LiquidFlow
}  // namespace ProcessLib
