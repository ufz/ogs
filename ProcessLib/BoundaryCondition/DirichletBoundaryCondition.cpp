/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "DirichletBoundaryCondition.h"

#include <algorithm>
#include <vector>
#include <logog/include/logog.hpp>

namespace ProcessLib
{
std::unique_ptr<DirichletBoundaryCondition> createDirichletBoundaryCondition(
    BaseLib::ConfigTree const& config, std::vector<std::size_t>&& mesh_node_ids,
    NumLib::LocalToGlobalIndexMap const& dof_table, std::size_t const mesh_id,
    int const variable_id, int const component_id)
{
    DBUG("Constructing DirichletBoundaryCondition from config.");
    //! \ogs_file_param{boundary_condition__type}
    config.checkConfigParameter("type", "Dirichlet");

    //! \ogs_file_param{boundary_condition__Dirichlet__value}
    auto const value = config.getConfigParameter<double>("value");
    DBUG("Using value %g", value);

    NumLib::IndexValueVector<GlobalIndexType> bc;

    // convert mesh node ids to global index for the given component
    bc.ids.reserve(bc.ids.size() + mesh_node_ids.size());
    bc.values.reserve(bc.values.size() + mesh_node_ids.size());
    for (auto& id : mesh_node_ids) {
        MeshLib::Location l(mesh_id, MeshLib::MeshItemType::Node, id);
        // TODO: that might be slow, but only done once
        const auto g_idx =
            dof_table.getGlobalIndex(l, variable_id, component_id);
        // For the DDC approach (e.g. with PETSc option), the negative
        // index of g_idx means that the entry by that index is a ghost one,
        // which should be dropped. Especially for PETSc routines MatZeroRows
        // and MatZeroRowsColumns, which are called to apply the Dirichlet BC,
        // the negative index is not accepted like other matrix or vector
        // PETSc routines. Therefore, the following if-condition is applied.
        if (g_idx >= 0) {
            bc.ids.emplace_back(g_idx);
            bc.values.emplace_back(value);
        }
    }

    return std::unique_ptr<DirichletBoundaryCondition>(
        new DirichletBoundaryCondition(std::move(bc)));
}

}  // namespace ProcessLib
