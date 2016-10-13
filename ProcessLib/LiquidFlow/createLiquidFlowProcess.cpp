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

    auto const gravity_param = config.getConfigParameter("gravitational_term");
    auto const axis =
        //! \ogs_file_attr{process__LIQUID_FLOW__gravitational_term__axis}
        gravity_param.getConfigAttributeOptional<std::string>("axis");
    // Gravitational acceleration
    auto const g =
        //! \ogs_file_attr{process__LIQUID_FLOW__gravitational_term__g}
        gravity_param.getConfigAttributeOptional<double>("g");

    int gravity_axis_id = -1;
    if (*axis == "x")
        gravity_axis_id = 0;
    else if (*axis == "y")
        gravity_axis_id = 1;
    else if (*axis == "z")
        gravity_axis_id = 2;

    //! \ogs_file_param{process__LIQUID_FLOW__material_property}
    auto const& mat_config = config.getConfigSubtree("material_property");

    auto const& mat_ids =
        mesh.getProperties().getPropertyVector<int>("MaterialIDs");
    if (mat_ids)
    {
        INFO("The liquid flow is in heterogeneous porous media.");
        return std::unique_ptr<Process>{new LiquidFlowProcess{
            mesh, std::move(jacobian_assembler), parameters, integration_order,
            std::move(process_variables), std::move(secondary_variables),
            std::move(named_function_caller), *mat_ids, gravity_axis_id, *g,
            mat_config}};
    }
    else
    {
        INFO("The liquid flow is in homogeneous porous media.");

        MeshLib::Properties dummy_property;
        auto const& dummy_property_vector =
            dummy_property.createNewPropertyVector<int>(
                "MaterialIDs", MeshLib::MeshItemType::Cell, 1);
        return std::unique_ptr<Process>{new LiquidFlowProcess{
            mesh, std::move(jacobian_assembler), parameters, integration_order,
            std::move(process_variables), std::move(secondary_variables),
            std::move(named_function_caller), *dummy_property_vector,
            gravity_axis_id, *g, mat_config}};
    }
}

}  // end of namespace
}  // end of namespace
