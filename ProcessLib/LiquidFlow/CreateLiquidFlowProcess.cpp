/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file   CreateLiquidFlowProcess.cpp
 *
 * Created on August 19, 2016, 1:30 PM
 */
#include "CreateLiquidFlowProcess.h"

#include <algorithm>

#include "MeshLib/MeshGenerators/MeshGenerator.h"

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
    //! \ogs_file_param{prj__processes__process__type}
    config.checkConfigParameter("type", "LIQUID_FLOW");

    DBUG("Create LiquidFlowProcess.");

    // Process variable.

    //! \ogs_file_param{prj__processes__process__LIQUID_FLOW__process_variables}
    auto const pv_config = config.getConfigSubtree("process_variables");

    auto process_variables = findProcessVariables(
        variables, pv_config,
        {//! \ogs_file_param_special{prj__processes__process__LIQUID_FLOW__process_variables__process_variable}
         "process_variable"});

    SecondaryVariableCollection secondary_variables;

    NumLib::NamedFunctionCaller named_function_caller({"LiquidFlow_pressure"});

    ProcessLib::parseSecondaryVariables(config, secondary_variables,
                                        named_function_caller);

    // Get the gravity vector for the Darcy velocity
    //! \ogs_file_param{prj__processes__process__LIQUID_FLOW__darcy_gravity}
    auto const& darcy_g_config = config.getConfigSubtree("darcy_gravity");
    const int gravity_axis_id_input =
        //! \ogs_file_param{prj__processes__process__LIQUID_FLOW__darcy_gravity__axis_id}
        darcy_g_config.getConfigParameter<int>("axis_id");
    assert(gravity_axis_id_input < static_cast<int>(mesh.getDimension()));
    const double g =
        //! \ogs_file_param{prj__processes__process__LIQUID_FLOW__darcy_gravity__g}
        darcy_g_config.getConfigParameter<double>("g");
    assert(g >= 0.);
    const int gravity_axis_id = (g == 0.) ? -1 : gravity_axis_id_input;

    //! \ogs_file_param{prj__processes__process__LIQUID_FLOW__material_property}
    auto const& mat_config = config.getConfigSubtree("material_property");

    auto const& mat_ids =
        mesh.getProperties().getPropertyVector<int>("MaterialIDs");
    if (mat_ids)
    {
        INFO("The liquid flow is in heterogeneous porous media.");
        const bool has_material_ids = true;
        return std::unique_ptr<Process>{new LiquidFlowProcess{
            mesh, std::move(jacobian_assembler), parameters, integration_order,
            std::move(process_variables), std::move(secondary_variables),
            std::move(named_function_caller), *mat_ids, has_material_ids,
            gravity_axis_id, g, mat_config}};
    }
    else
    {
        INFO("The liquid flow is in homogeneous porous media.");

        MeshLib::Properties dummy_property;
        // For a reference argument of LiquidFlowProcess(...).
        auto const& dummy_property_vector =
            dummy_property.createNewPropertyVector<int>(
                "MaterialIDs", MeshLib::MeshItemType::Cell, 1);

        // Since dummy_property_vector is only visible in this function,
        // the following constant, has_material_ids, is employed to indicate
        // that material_ids does not exist.
        const bool has_material_ids = false;

        return std::unique_ptr<Process>{new LiquidFlowProcess{
            mesh, std::move(jacobian_assembler), parameters, integration_order,
            std::move(process_variables), std::move(secondary_variables),
            std::move(named_function_caller), *dummy_property_vector,
            has_material_ids, gravity_axis_id, g, mat_config}};
    }
}

}  // end of namespace
}  // end of namespace
