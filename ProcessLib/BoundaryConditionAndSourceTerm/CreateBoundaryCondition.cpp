/**
 * \file
 * \copyright
 * Copyright (c) 2012-2023, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateBoundaryCondition.h"

#include "BaseLib/TimeInterval.h"
#include "BoundaryCondition.h"
#include "BoundaryConditionConfig.h"
#include "ConstraintDirichletBoundaryCondition.h"
#include "CreateDirichletBoundaryConditionWithinTimeInterval.h"
#include "DirichletBoundaryCondition.h"
#include "HCNonAdvectiveFreeComponentFlowBoundaryCondition.h"
#include "NeumannBoundaryCondition.h"
#include "NormalTractionBoundaryCondition.h"
#include "PhaseFieldIrreversibleDamageOracleBoundaryCondition.h"
#include "PrimaryVariableConstraintDirichletBoundaryCondition.h"
#include "Python/PythonBoundaryCondition.h"
#include "RobinBoundaryCondition.h"
#include "SolutionDependentDirichletBoundaryCondition.h"
#include "VariableDependentNeumannBoundaryCondition.h"

namespace ProcessLib
{
std::unique_ptr<BoundaryCondition> createBoundaryCondition(
    const BoundaryConditionConfig& config,
    const NumLib::LocalToGlobalIndexMap& dof_table,
    const MeshLib::Mesh& bulk_mesh, const int variable_id,
    const unsigned integration_order, const unsigned shapefunction_order,
    const std::vector<std::unique_ptr<ParameterLib::ParameterBase>>& parameters,
    const Process& process,
    [[maybe_unused]] std::vector<std::reference_wrapper<ProcessVariable>> const&
        all_process_variables_for_this_process,
    std::map<int,
             std::shared_ptr<MaterialPropertyLib::Medium>> const& /*media*/)
{
    // Surface mesh and bulk mesh must have equal axial symmetry flags!
    if (config.boundary_mesh.isAxiallySymmetric() !=
        bulk_mesh.isAxiallySymmetric())
    {
        OGS_FATAL(
            "The boundary mesh {:s} axially symmetric but the bulk mesh {:s}. "
            "Both must have an equal axial symmetry property.",
            config.boundary_mesh.isAxiallySymmetric() ? "is" : "is not",
            bulk_mesh.isAxiallySymmetric() ? "is" : "is not");
    }

    //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__type}
    auto const type = config.config.peekConfigParameter<std::string>("type");

    if (bool const component_id_required = type != "NormalTraction";
        component_id_required && !config.component_id.has_value())
    {
        OGS_FATAL(
            "Specifying the component id (<component>) for a boundary "
            "condition of type {} is mandatory.",
            type);
    }

    if (type == "Dirichlet")
    {
        return ProcessLib::createDirichletBoundaryCondition(
            config.config, config.boundary_mesh, dof_table, variable_id,
            *config.component_id, parameters);
    }
    if (type == "DirichletWithinTimeInterval")
    {
        return ProcessLib::createDirichletBoundaryConditionWithinTimeInterval(
            config.config, config.boundary_mesh, dof_table, variable_id,
            *config.component_id, parameters);
    }
    if (type == "Neumann")
    {
        return ProcessLib::createNeumannBoundaryCondition(
            config.config, config.boundary_mesh, dof_table, variable_id,
            *config.component_id, integration_order, shapefunction_order,
            bulk_mesh.getDimension(), parameters);
    }
    if (type == "Robin")
    {
        return ProcessLib::createRobinBoundaryCondition(
            config.config, config.boundary_mesh, dof_table, variable_id,
            *config.component_id, integration_order, shapefunction_order,
            bulk_mesh.getDimension(), parameters);
    }
    if (type == "VariableDependentNeumann")
    {
        return ProcessLib::createVariableDependentNeumannBoundaryCondition(
            config.config, config.boundary_mesh, dof_table, variable_id,
            *config.component_id, integration_order, shapefunction_order,
            bulk_mesh.getDimension(), parameters);
    }

    if (type == "Python")
    {
        return ProcessLib::createPythonBoundaryCondition(
            config.config, config.boundary_mesh, dof_table, bulk_mesh,
            variable_id, *config.component_id, integration_order,
            shapefunction_order, all_process_variables_for_this_process);
    }

    //
    // Special boundary conditions
    //
    if (type == "ConstraintDirichlet")
    {
        return createConstraintDirichletBoundaryCondition(
            config.config, config.boundary_mesh, dof_table, variable_id,
            integration_order, *config.component_id, parameters, process);
    }
    if (type == "PrimaryVariableConstraintDirichlet")
    {
        return createPrimaryVariableConstraintDirichletBoundaryCondition(
            config.config, config.boundary_mesh, dof_table, variable_id,
            *config.component_id, parameters);
    }
    if (type == "SolutionDependentDirichlet")
    {
        return ProcessLib::createSolutionDependentDirichletBoundaryCondition(
            config.config, config.boundary_mesh, dof_table, variable_id,
            *config.component_id, parameters);
    }
    if (type == "HCNonAdvectiveFreeComponentFlowBoundary")
    {
        return createHCNonAdvectiveFreeComponentFlowBoundaryCondition(
            config.config, config.boundary_mesh, dof_table, variable_id,
            *config.component_id, integration_order, parameters,
            bulk_mesh.getDimension(), process, shapefunction_order);
    }
    if (type == "NormalTraction")
    {
        switch (bulk_mesh.getDimension())
        {
            case 2:
                return ProcessLib::NormalTractionBoundaryCondition::
                    createNormalTractionBoundaryCondition<2>(
                        config.config, config.boundary_mesh, dof_table,
                        variable_id, integration_order, shapefunction_order,
                        parameters);
            case 3:
                return ProcessLib::NormalTractionBoundaryCondition::
                    createNormalTractionBoundaryCondition<3>(
                        config.config, config.boundary_mesh, dof_table,
                        variable_id, integration_order, shapefunction_order,
                        parameters);
            default:
                OGS_FATAL(
                    "NormalTractionBoundaryCondition can not be instantiated "
                    "for mesh dimensions other than two or three. "
                    "{}-dimensional mesh was given.",
                    bulk_mesh.getDimension());
        }
    }
    if (type == "PhaseFieldIrreversibleDamageOracleBoundaryCondition")
    {
        return ProcessLib::
            createPhaseFieldIrreversibleDamageOracleBoundaryCondition(
                //! \ogs_file_param_special{prj__process_variables__process_variable__boundary_conditions__boundary_condition__PhaseFieldIrreversibleDamageOracleBoundaryCondition}
                config.config, dof_table, bulk_mesh, variable_id,
                *config.component_id);
    }
    OGS_FATAL("Unknown boundary condition type: `{:s}'.", type);
}

}  // namespace ProcessLib
