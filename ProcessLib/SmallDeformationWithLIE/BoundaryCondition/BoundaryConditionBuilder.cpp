/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "BoundaryConditionBuilder.h"

#include "MeshGeoToolsLib/BoundaryElementsSearcher.h"
#include "MeshGeoToolsLib/MeshNodeSearcher.h"

#include "ProcessLib/BoundaryCondition/BoundaryConditionConfig.h"
#include "ProcessLib/BoundaryCondition/DirichletBoundaryCondition.h"
#include "ProcessLib/BoundaryCondition/RobinBoundaryCondition.h"

#include "NeumannBoundaryCondition.h"


static std::vector<MeshLib::Element*> getClonedElements(
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

namespace ProcessLib
{
namespace SmallDeformationWithLIE
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

        return createDirichletBoundaryCondition(
            config.config, std::move(ids), dof_table, mesh.getID(), variable_id,
            config.component_id, parameters);
    }
    else if (type == "Neumann")
    {
        return createNeumannBoundaryCondition(
            config.config,
            getClonedElements(boundary_element_searcher, config.geometry),
            dof_table, variable_id, config.component_id,
            mesh.isAxiallySymmetric(), integration_order, mesh.getDimension(),
            parameters, _fracture_prop);
    }
    else if (type == "Robin") {
        return createRobinBoundaryCondition(
            config.config,
            getClonedElements(boundary_element_searcher, config.geometry),
            dof_table, variable_id, config.component_id,
            mesh.isAxiallySymmetric(), integration_order, mesh.getDimension(),
            parameters);
    }
    else
    {
        OGS_FATAL("Unknown boundary condition type: `%s'.", type.c_str());
    }
}

}  // SmallDeformationWithLIE
}  // ProcessLib
