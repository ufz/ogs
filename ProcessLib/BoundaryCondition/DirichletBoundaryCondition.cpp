/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "DirichletBoundaryCondition.h"

#include <algorithm>
#include <logog/include/logog.hpp>
#include <vector>
#include "ProcessLib/Utils/ProcessUtils.h"

namespace ProcessLib
{
void DirichletBoundaryCondition::getEssentialBCValues(
    const double t, GlobalVector const& /*x*/,
    NumLib::IndexValueVector<GlobalIndexType>& bc_values) const
{
    SpatialPosition pos;

    auto const& bulk_node_ids_map =
        *_bc_mesh.getProperties().getPropertyVector<std::size_t>(
            "bulk_node_ids");

    bc_values.ids.clear();
    bc_values.values.clear();

    // convert mesh node ids to global index for the given component
    bc_values.ids.reserve(bc_values.ids.size() + _bc_mesh.getNumberOfNodes());
    bc_values.values.reserve(bc_values.values.size() +
                             _bc_mesh.getNumberOfNodes());
    for (auto const* const node : _bc_mesh.getNodes())
    {
        auto const id = bulk_node_ids_map[node->getID()];
        pos.setNodeID(id);
        MeshLib::Location l(_bulk_mesh_id, MeshLib::MeshItemType::Node, id);
        // TODO: that might be slow, but only done once
        const auto g_idx =
            _dof_table.getGlobalIndex(l, _variable_id, _component_id);
        if (g_idx == NumLib::MeshComponentMap::nop)
            continue;
        // For the DDC approach (e.g. with PETSc option), the negative
        // index of g_idx means that the entry by that index is a ghost one,
        // which should be dropped. Especially for PETSc routines MatZeroRows
        // and MatZeroRowsColumns, which are called to apply the Dirichlet BC,
        // the negative index is not accepted like other matrix or vector
        // PETSc routines. Therefore, the following if-condition is applied.
        if (g_idx >= 0)
        {
            bc_values.ids.emplace_back(g_idx);
            bc_values.values.emplace_back(_parameter(t, pos).front());
        }
    }
}

std::unique_ptr<DirichletBoundaryCondition> createDirichletBoundaryCondition(
    BaseLib::ConfigTree const& config, MeshLib::Mesh const& bc_mesh,
    NumLib::LocalToGlobalIndexMap const& dof_table,
    std::size_t const bulk_mesh_id, int const variable_id,
    int const component_id,
    const std::vector<std::unique_ptr<ProcessLib::ParameterBase>>& parameters)
{
    DBUG("Constructing DirichletBoundaryCondition from config.");
    //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__type}
    config.checkConfigParameter("type", "Dirichlet");

    //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__Dirichlet__parameter}
    auto const param_name = config.getConfigParameter<std::string>("parameter");
    DBUG("Using parameter %s", param_name.c_str());

    auto& param = findParameter<double>(param_name, parameters, 1);

    return std::make_unique<DirichletBoundaryCondition>(
        param, bc_mesh, dof_table, bulk_mesh_id, variable_id, component_id);
}

}  // namespace ProcessLib
