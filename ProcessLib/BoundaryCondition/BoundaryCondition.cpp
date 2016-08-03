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
#include "UniformNeumannBoundaryCondition.h"

namespace ProcessLib
{
std::unique_ptr<BoundaryCondition> createBoundaryCondition(
    const BoundaryConditionConfig& config,
    const NumLib::LocalToGlobalIndexMap& dof_table,
    const MeshLib::Mesh& mesh,
    const int variable_id,
    const unsigned integration_order)
{
    MeshGeoToolsLib::MeshNodeSearcher& mesh_node_searcher =
        MeshGeoToolsLib::MeshNodeSearcher::getMeshNodeSearcher(mesh);

    MeshGeoToolsLib::BoundaryElementsSearcher boundary_element_searcher(
        mesh, mesh_node_searcher);

    auto const type = config.config.peekConfigParameter<std::string>("type");

    if (type == "UniformDirichlet") {
        return createUniformDirichletBoundaryCondition(
            config, mesh_node_searcher, dof_table, variable_id);
    } else if (type == "UniformNeumann") {
        return createUniformNeumannBoundaryCondition(
            config, boundary_element_searcher, dof_table, variable_id,
            integration_order, mesh.getDimension());
    } else {
        OGS_FATAL("Unknown boundary condition type: `%s'.", type.c_str());
    }
}

}  // ProcessLib
