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
#include "MeshGeoToolsLib/BoundaryElementsSearcher.h"
#include "MeshGeoToolsLib/CreateSearchLength.h"
#include "MeshGeoToolsLib/MeshNodeSearcher.h"
#include "MeshGeoToolsLib/SearchLength.h"
#include "MeshLib/MeshEditing/DuplicateMeshComponents.h"
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
    std::unique_ptr<MeshGeoToolsLib::SearchLength> search_length_algorithm =
        MeshGeoToolsLib::createSearchLengthAlgorithm(config.config, mesh);

    MeshGeoToolsLib::MeshNodeSearcher const& mesh_node_searcher =
        MeshGeoToolsLib::MeshNodeSearcher::getMeshNodeSearcher(
            mesh, std::move(search_length_algorithm));

    // Find nodes' ids on the given mesh on which this boundary condition
    // is defined.
    std::vector<std::size_t> ids =
        mesh_node_searcher.getMeshNodeIDs(config.geometry);
    // Filter out ids, which are not part of mesh subsets corresponding to
    // the variable_id and component_id.

    // Sorted ids of all mesh_subsets.
    std::vector<std::size_t> sorted_nodes_ids;

    auto const& mesh_subset =
        dof_table.getMeshSubset(variable_id, *config.component_id);
    auto const& nodes = mesh_subset.getNodes();
    sorted_nodes_ids.reserve(sorted_nodes_ids.size() + nodes.size());
    std::transform(std::begin(nodes), std::end(nodes),
                   std::back_inserter(sorted_nodes_ids),
                   [](MeshLib::Node* const n) { return n->getID(); });
    std::sort(std::begin(sorted_nodes_ids), std::end(sorted_nodes_ids));

    auto ids_new_end_iterator = std::end(ids);
    ids_new_end_iterator = std::remove_if(
        std::begin(ids), ids_new_end_iterator,
        [&sorted_nodes_ids](std::size_t const node_id) {
            return !std::binary_search(std::begin(sorted_nodes_ids),
                                       std::end(sorted_nodes_ids), node_id);
        });
    ids.erase(ids_new_end_iterator, std::end(ids));

    DBUG("Found %d nodes for Dirichlet BCs for the variable %d and component %d",
         ids.size(), variable_id, *config.component_id);

    return ProcessLib::createDirichletBoundaryCondition(
        config.config, std::move(ids), dof_table, mesh.getID(), variable_id,
        *config.component_id, parameters);
}

std::unique_ptr<BoundaryCondition>
BoundaryConditionBuilder::createNeumannBoundaryCondition(
    const BoundaryConditionConfig& config,
    const NumLib::LocalToGlobalIndexMap& dof_table, const MeshLib::Mesh& mesh,
    const int variable_id, const unsigned integration_order,
    const unsigned shapefunction_order,
    const std::vector<std::unique_ptr<ProcessLib::ParameterBase>>& parameters)
{
    std::unique_ptr<MeshGeoToolsLib::SearchLength> search_length_algorithm =
        MeshGeoToolsLib::createSearchLengthAlgorithm(config.config, mesh);

    MeshGeoToolsLib::MeshNodeSearcher const& mesh_node_searcher =
        MeshGeoToolsLib::MeshNodeSearcher::getMeshNodeSearcher(
            mesh, std::move(search_length_algorithm));

    MeshGeoToolsLib::BoundaryElementsSearcher boundary_element_searcher(
        mesh, mesh_node_searcher);

    return ProcessLib::createNeumannBoundaryCondition(
        config.config,
        MeshLib::cloneElements(
            boundary_element_searcher.getBoundaryElements(config.geometry)),
        dof_table, variable_id, *config.component_id, mesh.isAxiallySymmetric(),
        integration_order, shapefunction_order, mesh.getDimension(),
        parameters);
}

std::unique_ptr<BoundaryCondition>
BoundaryConditionBuilder::createRobinBoundaryCondition(
    const BoundaryConditionConfig& config,
    const NumLib::LocalToGlobalIndexMap& dof_table, const MeshLib::Mesh& mesh,
    const int variable_id, const unsigned integration_order,
    const unsigned shapefunction_order,
    const std::vector<std::unique_ptr<ProcessLib::ParameterBase>>& parameters)
{
    std::unique_ptr<MeshGeoToolsLib::SearchLength> search_length_algorithm =
        MeshGeoToolsLib::createSearchLengthAlgorithm(config.config, mesh);

    MeshGeoToolsLib::MeshNodeSearcher const& mesh_node_searcher =
        MeshGeoToolsLib::MeshNodeSearcher::getMeshNodeSearcher(
            mesh, std::move(search_length_algorithm));

    MeshGeoToolsLib::BoundaryElementsSearcher boundary_element_searcher(
        mesh, mesh_node_searcher);

    return ProcessLib::createRobinBoundaryCondition(
        config.config,
        MeshLib::cloneElements(
            boundary_element_searcher.getBoundaryElements(config.geometry)),
        dof_table, variable_id, *config.component_id, mesh.isAxiallySymmetric(),
        integration_order, shapefunction_order, mesh.getDimension(),
        parameters);
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
    std::unique_ptr<MeshGeoToolsLib::SearchLength> search_length_algorithm =
        MeshGeoToolsLib::createSearchLengthAlgorithm(config.config, mesh);

    MeshGeoToolsLib::MeshNodeSearcher const& mesh_node_searcher =
        MeshGeoToolsLib::MeshNodeSearcher::getMeshNodeSearcher(
            mesh, std::move(search_length_algorithm));

    MeshGeoToolsLib::BoundaryElementsSearcher boundary_element_searcher(
        mesh, mesh_node_searcher);

    return ProcessLib::NormalTractionBoundaryCondition::
        createNormalTractionBoundaryCondition(
            config.config,
            MeshLib::cloneElements(
                boundary_element_searcher.getBoundaryElements(config.geometry)),
            dof_table, variable_id, mesh.isAxiallySymmetric(),
            integration_order, shapefunction_order, mesh.getDimension(),
            parameters);
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
