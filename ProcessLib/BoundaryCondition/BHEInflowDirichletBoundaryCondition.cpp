/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "BHEInflowDirichletBoundaryCondition.h"

#include <algorithm>
#include <logog/include/logog.hpp>
#include <vector>
#include "ProcessLib/Utils/ProcessUtils.h"

namespace ProcessLib
{
BHEInflowDirichletBoundaryCondition::BHEInflowDirichletBoundaryCondition(
    std::pair<GlobalIndexType, GlobalIndexType>&& in_out_global_indices,
    MeshLib::Mesh const& bc_mesh,
    std::vector<MeshLib::Node*> const& vec_inflow_bc_nodes,
    int const variable_id,
    int const component_id,
    std::unique_ptr<ProcessLib::HeatTransportBHE::BHE::BHEAbstract> const&
        pt_bhe)
    : _bc_mesh(bc_mesh), _pt_bhe(pt_bhe)
{
    DBUG(
        "Found %d nodes for BHE Inflow Dirichlet BCs for the variable %d and "
        "component %d",
        vec_inflow_bc_nodes.size(), variable_id, component_id);

    MeshLib::MeshSubset bc_mesh_subset{_bc_mesh, vec_inflow_bc_nodes};

    // create memory to store Tout values
    _T_out_values.clear();
    _T_out_indices.clear();

    _bc_values.ids.clear();
    _bc_values.values.clear();

    // convert mesh node ids to global index for the given component
    assert(bc_mesh_subset.getNumberOfNodes() == 1);
    _bc_values.ids.reserve(bc_mesh_subset.getNumberOfNodes());
    _bc_values.values.reserve(bc_mesh_subset.getNumberOfNodes());

    // that might be slow, but only done once
    const auto g_idx_T_in = in_out_global_indices.first;
    const auto g_idx_T_out = in_out_global_indices.second;

    if (g_idx_T_in >= 0 && g_idx_T_out >= 0)
    {
        _T_out_indices.emplace_back(g_idx_T_out);
        _T_out_values.emplace_back(320.0 /*using initial value*/);
        _bc_values.ids.emplace_back(g_idx_T_in);
        _bc_values.values.emplace_back(320.0 /*using initial value*/);
    }
}

void BHEInflowDirichletBoundaryCondition::getEssentialBCValues(
    const double t, GlobalVector const& x,
    NumLib::IndexValueVector<GlobalIndexType>& bc_values) const
{
    bc_values.ids.clear();
    bc_values.values.clear();

    bc_values.ids.resize(_bc_values.ids.size());
    bc_values.values.resize(_bc_values.values.size());

    const std::size_t n_nodes = _T_out_values.size();

    for (std::size_t i = 0; i < n_nodes; i++)
    {
        bc_values.ids[i] = _bc_values.ids[i];
        // here call the corresponding BHE functions
        auto const tmp_T_out = x[_T_out_indices[i]];
        bc_values.values[i] = _pt_bhe->getTinByTout(tmp_T_out, t);
    }
}

// update new values and corresponding indices.
void BHEInflowDirichletBoundaryCondition::preTimestep(const double /*t*/,
                                                      const GlobalVector& x)
{
    // for each BHE, the inflow temperature is dependent on
    // the ouflow temperature of the BHE.
    // Here the task is to get the outflow temperature and
    // save it locally
    auto const n_nodes = _bc_values.ids.size();
    for (std::size_t i = 0; i < n_nodes; i++)
    {
        // read the T_out
        _T_out_values[i] = x[_T_out_indices[i]];
    }
}

std::unique_ptr<BHEInflowDirichletBoundaryCondition>
createBHEInflowDirichletBoundaryCondition(
    std::pair<GlobalIndexType, GlobalIndexType>&& in_out_global_indices,
    MeshLib::Mesh const& bc_mesh,
    std::vector<MeshLib::Node*> const& vec_inflow_bc_nodes,
    int const variable_id, int const component_id,
    std::unique_ptr<ProcessLib::HeatTransportBHE::BHE::BHEAbstract> const&
        pt_bhe)
{
    DBUG("Constructing BHEInflowDirichletBoundaryCondition.");

    return std::make_unique<BHEInflowDirichletBoundaryCondition>(
        std::move(in_out_global_indices), bc_mesh, vec_inflow_bc_nodes,
        variable_id, component_id, pt_bhe);
}
}  // namespace ProcessLib