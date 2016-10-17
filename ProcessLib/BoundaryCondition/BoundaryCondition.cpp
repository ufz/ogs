/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "BoundaryCondition.h"
#include "MeshGeoToolsLib/BoundaryElementsSearcher.h"
#include "MeshGeoToolsLib/MeshNodeSearcher.h"
#include "BoundaryConditionConfig.h"
#include "DirichletBoundaryCondition.h"
#include "NeumannBoundaryCondition.h"
#include "RobinBoundaryCondition.h"

namespace ProcessLib
{

std::unique_ptr<BoundaryCondition> BoundaryConditionBuilder::createBoundaryCondition(
    const BoundaryConditionConfig& config,
    const NumLib::LocalToGlobalIndexMap& dof_table, const MeshLib::Mesh& mesh,
    const int variable_id, const unsigned integration_order,
    const std::vector<std::unique_ptr<ProcessLib::ParameterBase>>& parameters)
{
    MeshGeoToolsLib::MeshNodeSearcher& mesh_node_searcher =
        MeshGeoToolsLib::MeshNodeSearcher::getMeshNodeSearcher(mesh);

    MeshGeoToolsLib::BoundaryElementsSearcher boundary_element_searcher(
        mesh, mesh_node_searcher);

    //! \ogs_file_param{boundary_condition__type}
    auto const type = config.config.peekConfigParameter<std::string>("type");

    if (type == "Dirichlet")
    {
        return createDirichletBoundaryCondition(
                    config, dof_table, mesh, variable_id,
                    integration_order, parameters,
                    mesh_node_searcher, boundary_element_searcher);
    }
    else if (type == "Neumann")
    {
        return createNeumannBoundaryCondition(
                    config, dof_table, mesh, variable_id,
                    integration_order, parameters,
                    mesh_node_searcher, boundary_element_searcher);
    }
    else if (type == "Robin")
    {
        return createRobinBoundaryCondition(
                    config, dof_table, mesh, variable_id,
                    integration_order, parameters,
                    mesh_node_searcher, boundary_element_searcher);
    }
    else
    {
        OGS_FATAL("Unknown boundary condition type: `%s'.", type.c_str());
    }
}

std::unique_ptr<BoundaryCondition>
BoundaryConditionBuilder::createDirichletBoundaryCondition(
        const BoundaryConditionConfig& config,
        const NumLib::LocalToGlobalIndexMap& dof_table, const MeshLib::Mesh& mesh,
        const int variable_id, const unsigned /*integration_order*/,
        const std::vector<std::unique_ptr<ProcessLib::ParameterBase>>& parameters,
        MeshGeoToolsLib::MeshNodeSearcher& mesh_node_searcher,
        MeshGeoToolsLib::BoundaryElementsSearcher& /*boundary_element_searcher*/)
{
    // Find nodes' ids on the given mesh on which this boundary condition
    // is defined.
    std::vector<std::size_t> ids =
        mesh_node_searcher.getMeshNodeIDs(config.geometry);
    // Filter out ids, which are not part of mesh subsets corresponding to
    // the variable_id and component_id.

    // Sorted ids of all mesh_subsets.
    std::vector<std::size_t> sorted_nodes_ids;

    auto const& mesh_subsets =
        dof_table.getMeshSubsets(variable_id, config.component_id);
    for (auto const& mesh_subset : mesh_subsets)
    {
        auto const& nodes = mesh_subset->getNodes();
        sorted_nodes_ids.reserve(sorted_nodes_ids.size() + nodes.size());
        std::transform(std::begin(nodes), std::end(nodes),
                       std::back_inserter(sorted_nodes_ids),
                       [](MeshLib::Node* const n) { return n->getID(); });
    }
    std::sort(std::begin(sorted_nodes_ids), std::end(sorted_nodes_ids));

    auto ids_new_end_iterator = std::end(ids);
    ids_new_end_iterator = std::remove_if(
        std::begin(ids), ids_new_end_iterator,
        [&sorted_nodes_ids](std::size_t const node_id) {
            return !std::binary_search(std::begin(sorted_nodes_ids),
                                       std::end(sorted_nodes_ids), node_id);
        });
    ids.erase(ids_new_end_iterator, std::end(ids));

    return ProcessLib::createDirichletBoundaryCondition(
        config.config, std::move(ids), dof_table, mesh.getID(), variable_id,
        config.component_id, parameters);
}

std::unique_ptr<BoundaryCondition>
BoundaryConditionBuilder::createNeumannBoundaryCondition(
        const BoundaryConditionConfig& config,
        const NumLib::LocalToGlobalIndexMap& dof_table, const MeshLib::Mesh& mesh,
        const int variable_id, const unsigned integration_order,
        const std::vector<std::unique_ptr<ProcessLib::ParameterBase>>& parameters,
        MeshGeoToolsLib::MeshNodeSearcher& /*mesh_node_searcher*/,
        MeshGeoToolsLib::BoundaryElementsSearcher& boundary_element_searcher)
{
    return ProcessLib::createNeumannBoundaryCondition(
        config.config,
        getClonedElements(boundary_element_searcher, config.geometry),
        dof_table, variable_id, config.component_id,
        mesh.isAxiallySymmetric(), integration_order, mesh.getDimension(),
        parameters);
}

std::unique_ptr<BoundaryCondition>
BoundaryConditionBuilder::createRobinBoundaryCondition(
        const BoundaryConditionConfig& config,
        const NumLib::LocalToGlobalIndexMap& dof_table, const MeshLib::Mesh& mesh,
        const int variable_id, const unsigned integration_order,
        const std::vector<std::unique_ptr<ProcessLib::ParameterBase>>& parameters,
        MeshGeoToolsLib::MeshNodeSearcher& /*mesh_node_searcher*/,
        MeshGeoToolsLib::BoundaryElementsSearcher& boundary_element_searcher)
{
    return ProcessLib::createRobinBoundaryCondition(
        config.config,
        getClonedElements(boundary_element_searcher, config.geometry),
        dof_table, variable_id, config.component_id,
        mesh.isAxiallySymmetric(), integration_order, mesh.getDimension(),
        parameters);
}

std::vector<MeshLib::Element*> BoundaryConditionBuilder::getClonedElements(
    MeshGeoToolsLib::BoundaryElementsSearcher& boundary_element_searcher,
    GeoLib::GeoObject const& geometry)
{
    std::vector<MeshLib::Element*> elements =
        boundary_element_searcher.getBoundaryElements(geometry);

    // Deep copy all the elements, because the searcher might destroy the
    // originals. Store pointers to the copies in the elements vector (i.e.,
    // in-place modification).
    for (auto& e : elements)
        e = e->clone();

    return elements;
}

}  // ProcessLib
