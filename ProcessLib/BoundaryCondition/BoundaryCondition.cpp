/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "BoundaryCondition.h"
#include "BoundaryConditionConfig.h"
#include "DirichletBoundaryCondition.h"
#include "ConstraintDirichletBoundaryCondition.h"
#include "NeumannBoundaryCondition.h"
#include "NonuniformDirichletBoundaryCondition.h"
#include "NonuniformNeumannBoundaryCondition.h"
#include "NormalTractionBoundaryCondition.h"
#include "PhaseFieldIrreversibleDamageOracleBoundaryCondition.h"
#include "RobinBoundaryCondition.h"

namespace ProcessLib
{
std::unique_ptr<BoundaryCondition>
BoundaryConditionBuilder::createBoundaryCondition(
    const BoundaryConditionConfig& config,
    const NumLib::LocalToGlobalIndexMap& dof_table, const MeshLib::Mesh& mesh,
    const int variable_id, const unsigned integration_order,
    const unsigned shapefunction_order,
    const std::vector<std::unique_ptr<ProcessLib::ParameterBase>>& parameters)
{
    //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__type}
    auto const type = config.config.peekConfigParameter<std::string>("type");

    if (type == "Dirichlet")
    {
        return createDirichletBoundaryCondition(
                    config, dof_table, mesh, variable_id,
                    integration_order, shapefunction_order, parameters);
    }
    if (type == "Neumann")
    {
        return createNeumannBoundaryCondition(
                    config, dof_table, mesh, variable_id,
                    integration_order, shapefunction_order, parameters);
    }
    if (type == "Robin")
    {
        return createRobinBoundaryCondition(
                    config, dof_table, mesh, variable_id,
                    integration_order, shapefunction_order, parameters);
    }
    if (type == "NonuniformDirichlet")
    {
        return createNonuniformDirichletBoundaryCondition(config, dof_table,
                                                          mesh, variable_id);
    }
    if (type == "NonuniformNeumann")
    {
        return createNonuniformNeumannBoundaryCondition(
            config, dof_table, mesh, variable_id, integration_order,
            shapefunction_order);
    }
    //
    // Special boundary conditions
    //
    if (type == "NormalTraction")
    {
        return createNormalTractionBoundaryCondition(
            config, dof_table, mesh, variable_id, integration_order,
            shapefunction_order, parameters);
    }
    if (type == "PhaseFieldIrreversibleDamageOracleBoundaryCondition")
    {
        return createPhaseFieldIrreversibleDamageOracleBoundaryCondition(
            config, dof_table, mesh, variable_id, integration_order,
            shapefunction_order, parameters);
    }
    OGS_FATAL("Unknown boundary condition type: `%s'.", type.c_str());
}

std::unique_ptr<BoundaryCondition>
BoundaryConditionBuilder::createDirichletBoundaryCondition(
    const BoundaryConditionConfig& config,
    const NumLib::LocalToGlobalIndexMap& dof_table, const MeshLib::Mesh& mesh,
    const int variable_id, const unsigned /*integration_order*/,
    const unsigned /*shapefunction_order*/,
    const std::vector<std::unique_ptr<ProcessLib::ParameterBase>>& parameters)
{
    return ProcessLib::createDirichletBoundaryCondition(
        config.config, config.mesh, dof_table, mesh.getID(), variable_id,
        *config.component_id, parameters);
}

std::unique_ptr<BoundaryCondition>
BoundaryConditionBuilder::createConstraintDirichletBoundaryCondition(
    const BoundaryConditionConfig& config,
    const NumLib::LocalToGlobalIndexMap& dof_table, const MeshLib::Mesh& mesh,
    const int variable_id, const unsigned integration_order,
    const std::vector<std::unique_ptr<ProcessLib::ParameterBase>>& parameters,
    Process const& process)
{
    return ProcessLib::createConstraintDirichletBoundaryCondition(
        config.config, config.mesh, dof_table, variable_id, integration_order,
        *config.component_id, parameters, process);
}

std::unique_ptr<BoundaryCondition>
BoundaryConditionBuilder::createNeumannBoundaryCondition(
    const BoundaryConditionConfig& config,
    const NumLib::LocalToGlobalIndexMap& dof_table, const MeshLib::Mesh& mesh,
    const int variable_id, const unsigned integration_order,
    const unsigned shapefunction_order,
    const std::vector<std::unique_ptr<ProcessLib::ParameterBase>>& parameters)
{
    return ProcessLib::createNeumannBoundaryCondition(
        config.config, config.mesh, dof_table, variable_id,
        *config.component_id, mesh.isAxiallySymmetric(), integration_order,
        shapefunction_order, mesh.getDimension(), parameters);
}

std::unique_ptr<BoundaryCondition>
BoundaryConditionBuilder::createRobinBoundaryCondition(
    const BoundaryConditionConfig& config,
    const NumLib::LocalToGlobalIndexMap& dof_table, const MeshLib::Mesh& mesh,
    const int variable_id, const unsigned integration_order,
    const unsigned shapefunction_order,
    const std::vector<std::unique_ptr<ProcessLib::ParameterBase>>& parameters)
{
    return ProcessLib::createRobinBoundaryCondition(
        config.config, config.mesh, dof_table, variable_id,
        *config.component_id, mesh.isAxiallySymmetric(), integration_order,
        shapefunction_order, mesh.getDimension(), parameters);
}

std::unique_ptr<BoundaryCondition>
BoundaryConditionBuilder::createNonuniformDirichletBoundaryCondition(
    const BoundaryConditionConfig& config,
    const NumLib::LocalToGlobalIndexMap& dof_table, const MeshLib::Mesh& mesh,
    const int variable_id)
{
    return ProcessLib::createNonuniformDirichletBoundaryCondition(
        config.config, dof_table, variable_id, *config.component_id, mesh);
}

std::unique_ptr<BoundaryCondition>
BoundaryConditionBuilder::createNonuniformNeumannBoundaryCondition(
    const BoundaryConditionConfig& config,
    const NumLib::LocalToGlobalIndexMap& dof_table, const MeshLib::Mesh& mesh,
    const int variable_id, const unsigned integration_order,
    const unsigned shapefunction_order)
{
    return ProcessLib::createNonuniformNeumannBoundaryCondition(
        config.config, dof_table, variable_id, *config.component_id,
        integration_order, shapefunction_order, mesh);
}

std::unique_ptr<BoundaryCondition>
BoundaryConditionBuilder::createNormalTractionBoundaryCondition(
    const BoundaryConditionConfig& config,
    const NumLib::LocalToGlobalIndexMap& dof_table, const MeshLib::Mesh& mesh,
    const int variable_id, const unsigned integration_order,
    const unsigned shapefunction_order,
    const std::vector<std::unique_ptr<ProcessLib::ParameterBase>>& parameters)
{
    return ProcessLib::NormalTractionBoundaryCondition::
        createNormalTractionBoundaryCondition(
            config.config, config.mesh, dof_table, variable_id,
            mesh.isAxiallySymmetric(), integration_order, shapefunction_order,
            mesh.getDimension(), parameters);
}

std::unique_ptr<BoundaryCondition> BoundaryConditionBuilder::
    createPhaseFieldIrreversibleDamageOracleBoundaryCondition(
        const BoundaryConditionConfig& config,
        const NumLib::LocalToGlobalIndexMap& dof_table,
        const MeshLib::Mesh& mesh, const int variable_id,
        const unsigned /*integration_order*/,
        const unsigned /*shapefunction_order*/,
        const std::vector<
            std::unique_ptr<ProcessLib::ParameterBase>>& /*parameters*/)
{
    return ProcessLib::
        createPhaseFieldIrreversibleDamageOracleBoundaryCondition(
            config.config, dof_table, mesh, variable_id, *config.component_id);
}

}  // namespace ProcessLib
