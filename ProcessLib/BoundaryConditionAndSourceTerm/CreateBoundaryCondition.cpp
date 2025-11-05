/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateBoundaryCondition.h"

#include "BaseLib/Error.h"
#include "BaseLib/TimeInterval.h"
#include "BoundaryCondition.h"
#include "BoundaryConditionConfig.h"
#include "ConstraintDirichletBoundaryCondition.h"
#include "CreateDirichletBoundaryConditionWithinTimeInterval.h"
#include "CreateReleaseNodalForce.h"
#include "CreateTimeDecayDirichletBoundaryCondition.h"
#include "DirichletBoundaryCondition.h"
#include "HCNonAdvectiveFreeComponentFlowBoundaryCondition.h"
#include "NeumannBoundaryCondition.h"
#include "NormalTractionBoundaryCondition.h"
#include "ParameterLib/Parameter.h"
#include "PhaseFieldIrreversibleDamageOracleBoundaryCondition.h"
#include "PrimaryVariableConstraintDirichletBoundaryCondition.h"
#include "Python/PythonBoundaryCondition.h"
#include "RobinBoundaryCondition.h"
#include "SolutionDependentDirichletBoundaryCondition.h"
#include "VariableDependentNeumannBoundaryCondition.h"
#include "WellboreCompensateNeumannBoundaryCondition.h"

namespace ProcessLib
{
std::vector<std::unique_ptr<BoundaryCondition>> createBoundaryCondition(
    const BoundaryConditionConfig& config,
    const NumLib::LocalToGlobalIndexMap& dof_table,
    const MeshLib::Mesh& bulk_mesh, const int variable_id,
    const unsigned integration_order, const unsigned shapefunction_order,
    const std::vector<std::unique_ptr<ParameterLib::ParameterBase>>& parameters,
    const Process& process,
    [[maybe_unused]] std::vector<std::reference_wrapper<ProcessVariable>> const&
        all_process_variables_for_this_process,
    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const& media)
{
    // Surface mesh and bulk mesh must have equal axial symmetry flags!
    for (auto const& i : config.boundary_meshes)
    {
        if (i.get().isAxiallySymmetric() != bulk_mesh.isAxiallySymmetric())
        {
            OGS_FATAL(
                "The boundary mesh {:s} axially symmetric but the bulk mesh "
                "{:s}. "
                "Both must have an equal axial symmetry property.",
                i.get().isAxiallySymmetric() ? "is" : "is not",
                bulk_mesh.isAxiallySymmetric() ? "is" : "is not");
        }
    }

    std::vector<std::unique_ptr<BoundaryCondition>> boundary_conditions;

    //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__type}
    auto const type = config.config.peekConfigParameter<std::string>("type");

    if (bool const component_id_required =
            !(type == "NormalTraction" || type == "ReleaseNodalForce");
        component_id_required && !config.component_id.has_value())
    {
        OGS_FATAL(
            "Specifying the component id (<component>) for a boundary "
            "condition of type {} is mandatory.",
            type);
    }

    if (type == "Dirichlet")
    {
        auto const args = parseDirichletBCConfig(config.config);
        std::vector<std::unique_ptr<BoundaryCondition>> dirichlet_conditions;
        for (auto const& bc_mesh : config.boundary_meshes)
        {
            dirichlet_conditions.push_back(
                ProcessLib::createDirichletBoundaryCondition(
                    args, bc_mesh, dof_table, variable_id, *config.component_id,
                    parameters));
        }
        return dirichlet_conditions;
    }
    if (type == "DirichletWithinTimeInterval")
    {
        auto const args =
            parseDirichletBoundaryConditionWithinTimeIntervalConfig(
                config.config);
        std::vector<std::unique_ptr<BoundaryCondition>> conditions;
        for (auto const& bc_mesh : config.boundary_meshes)
        {
            conditions.push_back(
                ProcessLib::createDirichletBoundaryConditionWithinTimeInterval(
                    args.first, args.second, bc_mesh, dof_table, variable_id,
                    *config.component_id, parameters));
        }
        return conditions;
    }
    if (type == "TimeDecayDirichlet")
    {
        auto const [parameter_name, lower_limit] =
            parseTimeDecayDirichletBoundaryConditionConfig(config.config);
        std::vector<std::unique_ptr<BoundaryCondition>> conditions;
        for (auto const& bc_mesh : config.boundary_meshes)
        {
            conditions.push_back(
                ProcessLib::createTimeDecayDirichletBoundaryCondition(
                    variable_id, *config.component_id, parameter_name,
                    lower_limit, bc_mesh, dof_table, parameters));
        }
        return conditions;
    }
    if (type == "Neumann")
    {
        auto const [parameter_name, area_parameter_name] =
            parseNeumannBoundaryCondition(config.config);
        std::vector<std::unique_ptr<BoundaryCondition>> conditions;
        for (auto const& bc_mesh : config.boundary_meshes)
        {
            conditions.push_back(ProcessLib::createNeumannBoundaryCondition(
                parameter_name, area_parameter_name, bc_mesh, dof_table,
                variable_id, *config.component_id, integration_order,
                shapefunction_order, bulk_mesh.getDimension(), parameters));
        }
        return conditions;
    }
    if (type == "Robin")
    {
        auto const [alpha_name, u_0_name, area_parameter_name] =
            parseRobinBoundaryCondition(config.config);
        std::vector<std::unique_ptr<BoundaryCondition>> conditions;
        for (auto const& bc_mesh : config.boundary_meshes)
        {
            conditions.push_back(ProcessLib::createRobinBoundaryCondition(
                alpha_name, u_0_name, area_parameter_name, bc_mesh, dof_table,
                variable_id, *config.component_id, integration_order,
                shapefunction_order, bulk_mesh.getDimension(), parameters));
        }
        return conditions;
    }
    if (type == "VariableDependentNeumann")
    {
        auto const [constant_name, coefficient_current_variable_name,
                    coefficient_other_variable_name,
                    coefficient_mixed_variables_name] =
            parseVariableDependentNeumannBoundaryCondition(config.config);
        std::vector<std::unique_ptr<BoundaryCondition>> conditions;
        for (auto const& bc_mesh : config.boundary_meshes)
        {
            conditions.push_back(
                ProcessLib::createVariableDependentNeumannBoundaryCondition(
                    constant_name, coefficient_current_variable_name,
                    coefficient_other_variable_name,
                    coefficient_mixed_variables_name, bc_mesh, dof_table,
                    variable_id, *config.component_id, integration_order,
                    shapefunction_order, bulk_mesh.getDimension(), parameters));
        }
        return conditions;
    }
    if (type == "Python")
    {
        auto const [bc_object, flush_stdout] =
            parsePythonBoundaryCondition(config.config);
        std::vector<std::unique_ptr<BoundaryCondition>> conditions;
        for (auto const& bc_mesh : config.boundary_meshes)
        {
            conditions.push_back(ProcessLib::createPythonBoundaryCondition(
                bc_object, flush_stdout, bc_mesh, dof_table, bulk_mesh,
                variable_id, *config.component_id, integration_order,
                shapefunction_order, all_process_variables_for_this_process));
        }
        return conditions;
    }

    //
    // Special boundary conditions
    //
    if (type == "ConstraintDirichlet")
    {
        auto const [type, process_variable, threshold, direction_string, name] =
            parseConstraintDirichletBoundaryCondition(config.config);
        std::vector<std::unique_ptr<BoundaryCondition>> conditions;
        for (auto const& bc_mesh : config.boundary_meshes)
        {
            conditions.push_back(
                ProcessLib::createConstraintDirichletBoundaryCondition(
                    type, process_variable, threshold, direction_string, name,
                    bc_mesh, dof_table, variable_id, integration_order,
                    *config.component_id, parameters, process));
        }
        return conditions;
    }
    if (type == "PrimaryVariableConstraintDirichlet")
    {
        auto const [name, threshold_parameter_name, comparison_operator] =
            parsePrimaryVariableConstraintDirichletBoundaryCondition(
                config.config);
        std::vector<std::unique_ptr<BoundaryCondition>> conditions;
        for (auto const& bc_mesh : config.boundary_meshes)
        {
            conditions.push_back(
                createPrimaryVariableConstraintDirichletBoundaryCondition(
                    name, threshold_parameter_name, comparison_operator,
                    bc_mesh, dof_table, variable_id, *config.component_id,
                    parameters));
        }
        return conditions;
    }
    if (type == "WellboreCompensateNeumann")
    {
        auto const coefficients =
            parseWellboreCompensateNeumannBoundaryCondition(config.config);
        std::vector<std::unique_ptr<BoundaryCondition>> conditions;
        for (auto const& bc_mesh : config.boundary_meshes)
        {
            conditions.push_back(
                ProcessLib::createWellboreCompensateNeumannBoundaryCondition(
                    coefficients, bc_mesh, dof_table, variable_id,
                    *config.component_id, integration_order,
                    shapefunction_order, bulk_mesh.getDimension(), media));
        }
        return conditions;
    }
    if (type == "SolutionDependentDirichlet")
    {
        auto const [property_name, initial_value_parameter_string] =
            parseSolutionDependentDirichletBoundaryCondition(config.config);
        std::vector<std::unique_ptr<BoundaryCondition>> conditions;
        for (auto const& bc_mesh : config.boundary_meshes)
        {
            conditions.push_back(
                ProcessLib::createSolutionDependentDirichletBoundaryCondition(
                    property_name, initial_value_parameter_string, bc_mesh,
                    dof_table, variable_id, *config.component_id, parameters));
        }
        return conditions;
    }
    if (type == "HCNonAdvectiveFreeComponentFlowBoundary")
    {
        auto const boundary_permeability_name =
            parseHCNonAdvectiveFreeComponentFlowBoundaryCondition(
                config.config);
        std::vector<std::unique_ptr<BoundaryCondition>> conditions;
        for (auto const& bc_mesh : config.boundary_meshes)
        {
            conditions.push_back(
                createHCNonAdvectiveFreeComponentFlowBoundaryCondition(
                    boundary_permeability_name, bc_mesh, dof_table, variable_id,
                    *config.component_id, integration_order, parameters,
                    bulk_mesh.getDimension(), process, shapefunction_order));
        }
        return conditions;
    }
    if (type == "NormalTraction")
    {
        auto const parameter_name = NormalTractionBoundaryCondition::
            parseNormalTractionBoundaryCondition(config.config);
        std::vector<std::unique_ptr<BoundaryCondition>> conditions;
        for (auto const& bc_mesh : config.boundary_meshes)
        {
            switch (bulk_mesh.getDimension())
            {
                case 2:
                    conditions.push_back(
                        ProcessLib::NormalTractionBoundaryCondition::
                            createNormalTractionBoundaryCondition<2>(
                                parameter_name, bc_mesh, bulk_mesh, dof_table,
                                variable_id, integration_order,
                                shapefunction_order, parameters));
                    break;
                case 3:
                    conditions.push_back(
                        ProcessLib::NormalTractionBoundaryCondition::
                            createNormalTractionBoundaryCondition<3>(
                                parameter_name, bc_mesh, bulk_mesh, dof_table,
                                variable_id, integration_order,
                                shapefunction_order, parameters));
                    break;
                default:
                    OGS_FATAL(
                        "NormalTractionBoundaryCondition can not be "
                        "instantiated for mesh dimensions other than two or "
                        "three. {}-dimensional mesh was given.",
                        bulk_mesh.getDimension());
            }
        }
        return conditions;
    }
    if (type == "PhaseFieldIrreversibleDamageOracleBoundaryCondition")
    {
        parsePhaseFieldIrreversibleDamageOracleBoundaryCondition(config.config);
        std::vector<std::unique_ptr<BoundaryCondition>> conditions;
        conditions.push_back(
            createPhaseFieldIrreversibleDamageOracleBoundaryCondition(
                //! \ogs_file_param_special{prj__process_variables__process_variable__boundary_conditions__boundary_condition__PhaseFieldIrreversibleDamageOracleBoundaryCondition}
                dof_table, bulk_mesh, variable_id, *config.component_id));
        return conditions;
    }
    if (type == "ReleaseNodalForce")
    {
        auto const decay_parameter_name = parseReleaseNodalForce(config);
        std::vector<std::unique_ptr<BoundaryCondition>> conditions;
        for (auto const& bc_mesh : config.boundary_meshes)
        {
            conditions.push_back(ProcessLib::createReleaseNodalForce(
                //! \ogs_file_param_special{prj__process_variables__process_variable__boundary_conditions__boundary_condition__ReleaseNodalForce}
                bulk_mesh.getDimension(), variable_id, decay_parameter_name,
                config, bc_mesh, dof_table, parameters));
        }
        return conditions;
    }
    OGS_FATAL("Unknown boundary condition type: `{:s}'.", type);
}

}  // namespace ProcessLib
