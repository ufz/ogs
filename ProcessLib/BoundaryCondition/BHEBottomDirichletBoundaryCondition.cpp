/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "BHEBottomDirichletBoundaryCondition.h"

#include <algorithm>
#include <logog/include/logog.hpp>
#include <vector>
#include "ProcessLib/Utils/ProcessUtils.h"

namespace ProcessLib
{
BHEBottomDirichletBoundaryCondition::BHEBottomDirichletBoundaryCondition(
    GlobalIndexType const global_idx_T_in_bottom,
    GlobalIndexType const global_idx_T_out_bottom,
    MeshLib::Mesh const& bulk_mesh,
    std::vector<MeshLib::Node*> const& vec_outflow_bc_nodes,
    int const variable_id,
    int const component_id)
    : _bulk_mesh(bulk_mesh)
{
    DBUG(
        "Found %d nodes for BHE bottom Dirichlet BCs for the variable %d "
        "and "
        "component %d",
        vec_outflow_bc_nodes.size(), variable_id, component_id);

    MeshLib::MeshSubset bc_mesh_subset{_bulk_mesh, vec_outflow_bc_nodes};

    // create memory to store Tout values
    _T_in_values.ids.clear();
    _T_in_values.values.clear();

    _bc_values.ids.clear();
    _bc_values.values.clear();

    // convert mesh node ids to global index for the given component
    _bc_values.ids.reserve(bc_mesh_subset.getNumberOfNodes());
    _bc_values.values.reserve(bc_mesh_subset.getNumberOfNodes());

    const auto g_T_out_idx = global_idx_T_out_bottom;
    if (g_T_out_idx >= 0)
    {
        _bc_values.ids.emplace_back(g_T_out_idx);
        _bc_values.values.emplace_back(298.15);
    }

    const auto g_T_in_idx = global_idx_T_in_bottom;

    if (g_T_in_idx >= 0)
    {
        _T_in_values.ids.emplace_back(g_T_in_idx);
        _T_in_values.values.emplace_back(298.15);
    }
}

void BHEBottomDirichletBoundaryCondition::getEssentialBCValues(
    const double /*t*/, GlobalVector const& x,
    NumLib::IndexValueVector<GlobalIndexType>& bc_values) const
{
    bc_values.ids.clear();
    bc_values.values.clear();

    bc_values.ids.resize(_bc_values.ids.size());
    bc_values.values.resize(_bc_values.values.size());

    const std::size_t n_nodes = _T_in_values.ids.size();
    for (std::size_t i = 0; i < n_nodes; i++)
    {
        bc_values.ids[i] = _bc_values.ids[i];
        // here, the outflow temperature is always
        // the same as the inflow temperature
        // get the inflow temperature from here.
        bc_values.values[i] = x[_T_in_values.ids[i]];
    }
}

// update new values and corresponding indices.
void BHEBottomDirichletBoundaryCondition::preTimestep(const double /*t*/,
                                                      const GlobalVector& x)
{
    // At the bottom of each BHE, the outflow temperature
    // is the same as the inflow temperature.
    // Here the task is to get the inflow temperature and
    // save it locally
    auto const n_nodes = _bc_values.ids.size();
    for (std::size_t i = 0; i < n_nodes; i++)
    {
        // read the T_out
        _T_in_values.values[i] = x[_T_in_values.ids[i]];
    }
}

std::unique_ptr<BHEBottomDirichletBoundaryCondition>
createBHEBottomDirichletBoundaryCondition(
    GlobalIndexType global_idx_T_in_bottom,
    GlobalIndexType global_idx_T_out_bottom, MeshLib::Mesh const& bulk_mesh,
    std::vector<MeshLib::Node*> const& vec_outflow_bc_nodes,
    int const variable_id, int const component_id)
{
    DBUG("Constructing BHEBottomDirichletBoundaryCondition from config.");

    return std::make_unique<BHEBottomDirichletBoundaryCondition>(
        global_idx_T_in_bottom, global_idx_T_out_bottom, bulk_mesh,
        vec_outflow_bc_nodes, variable_id, component_id);
}
}  // namespace ProcessLib
