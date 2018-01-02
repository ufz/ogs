/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "SourceTermBuilder.h"
#include "SourceTermConfig.h"
#include "CreateNodalSourceTerm.h"
#include "NodalSourceTerm.h"
#include "MeshGeoToolsLib/BoundaryElementsSearcher.h"
#include "MeshGeoToolsLib/CreateSearchLength.h"
#include "MeshGeoToolsLib/MeshNodeSearcher.h"
#include "MeshGeoToolsLib/SearchLength.h"

namespace ProcessLib
{
std::unique_ptr<NodalSourceTerm> SourceTermBuilder::createSourceTerm(
    const SourceTermConfig& config,
    const NumLib::LocalToGlobalIndexMap& dof_table, const MeshLib::Mesh& mesh,
    const int variable_id, const unsigned integration_order,
    const unsigned shapefunction_order)
{
    //! \ogs_file_param{prj__process_variables__process_variable__source_terms__source_term__type}
    auto const type = config.config.peekConfigParameter<std::string>("type");

    if (type == "Nodal")
    {
        return createNodalSourceTerm(config, dof_table, mesh, variable_id,
                                     integration_order, shapefunction_order);
    }

    OGS_FATAL("Unknown source term type: `%s'.", type.c_str());
}

std::unique_ptr<NodalSourceTerm> SourceTermBuilder::createNodalSourceTerm(
    const SourceTermConfig& config,
    const NumLib::LocalToGlobalIndexMap& dof_table, const MeshLib::Mesh& mesh,
    const int variable_id, const unsigned /*integration_order*/,
    const unsigned /*shapefunction_order*/)
{
    std::unique_ptr<MeshGeoToolsLib::SearchLength> search_length_algorithm =
        MeshGeoToolsLib::createSearchLengthAlgorithm(config.config, mesh);

    MeshGeoToolsLib::MeshNodeSearcher const& mesh_node_searcher =
        MeshGeoToolsLib::MeshNodeSearcher::getMeshNodeSearcher(
            mesh, std::move(search_length_algorithm));

    // Find nodes' ids on the given mesh on which this source term is defined.
    std::vector<std::size_t> ids =
        mesh_node_searcher.getMeshNodeIDs(config.geometry);

    // Filter out ids, which are not part of mesh subsets corresponding to
    // the variable_id and component_id.

    // Sorted ids of all mesh_subsets.
    std::vector<std::size_t> sorted_nodes_ids;

    auto const& mesh_subsets =
        dof_table.getMeshSubsets(variable_id, *config.component_id);
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

    DBUG(
        "Found %d nodes for nodal source term for the variable %d and "
        "component %d",
        ids.size(), variable_id, *config.component_id);

    if (ids.size() != 1)
        OGS_FATAL(
            "Found %d nodes for nodal source term, but exactly one node is "
            "required.");

    return ProcessLib::createNodalSourceTerm(
        config.config, dof_table, mesh.getID(), ids[0], variable_id,
        *config.component_id);
}

}  // ProcessLib
