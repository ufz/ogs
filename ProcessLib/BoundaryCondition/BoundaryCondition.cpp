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
#include "UniformDirichletBoundaryCondition.h"
#include "NeumannBoundaryCondition.h"
#include "RobinBoundaryCondition.h"

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
std::unique_ptr<BoundaryCondition> createBoundaryCondition(
    const BoundaryConditionConfig& config,
    const NumLib::LocalToGlobalIndexMap& dof_table,
    const MeshLib::Mesh& mesh,
    const int variable_id,
    const unsigned integration_order,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters)
{
    MeshGeoToolsLib::MeshNodeSearcher& mesh_node_searcher =
        MeshGeoToolsLib::MeshNodeSearcher::getMeshNodeSearcher(mesh);

    MeshGeoToolsLib::BoundaryElementsSearcher boundary_element_searcher(
        mesh, mesh_node_searcher);

    //! \ogs_file_param{boundary_condition__type}
    auto const type = config.config.peekConfigParameter<std::string>("type");

    if (type == "UniformDirichlet")
    {
        // Find nodes' ids on the given mesh on which this boundary condition
        // is defined.
        std::vector<std::size_t> ids =
            mesh_node_searcher.getMeshNodeIDs(config.geometry);

        return createUniformDirichletBoundaryCondition(
            config.config, std::move(ids), dof_table, mesh.getID(), variable_id,
            config.component_id);
    }
    else if (type == "Neumann")
    {
        return createNeumannBoundaryCondition(
            config.config,
            getClonedElements(boundary_element_searcher, config.geometry),
            dof_table, variable_id, config.component_id, integration_order,
            mesh.getDimension(), parameters);
    }
    else if (type == "Robin") {
        return createRobinBoundaryCondition(
            config.config,
            getClonedElements(boundary_element_searcher, config.geometry),
            dof_table, variable_id, config.component_id, integration_order,
            mesh.getDimension(), parameters);
    }
    else
    {
        OGS_FATAL("Unknown boundary condition type: `%s'.", type.c_str());
    }
}

}  // ProcessLib
