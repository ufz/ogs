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
#include "ProcessLib/Utils/ProcessUtils.h"

namespace ProcessLib
{
void DirichletBoundaryCondition::preTimestep(const double /*t*/)
{
    if (_parameter.isTimeDependent())
        _already_computed = false;
}

void DirichletBoundaryCondition::getEssentialBCValues(
    const double t, NumLib::IndexValueVector<GlobalIndexType>& bc_values) const
{
    // TODO: Reenable when fixed ;)
    //if (_already_computed)
        //return;

    _already_computed = true;

    SpatialPosition pos;

    bc_values.ids.clear();
    bc_values.values.clear();

    // convert mesh node ids to global index for the given component
    bc_values.ids.reserve(bc_values.ids.size() + _mesh_node_ids.size());
    bc_values.values.reserve(bc_values.values.size() + _mesh_node_ids.size());
    for (auto const id : _mesh_node_ids)
    {
        pos.setNodeID(id);
        MeshLib::Location l(_mesh_id, MeshLib::MeshItemType::Node, id);
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
        if (g_idx >= 0) {
            bc_values.ids.emplace_back(g_idx);
            bc_values.values.emplace_back(_parameter(t, pos).front());
        }
    }
}

std::unique_ptr<DirichletBoundaryCondition> createDirichletBoundaryCondition(
    BaseLib::ConfigTree const& config, std::vector<std::size_t>&& mesh_node_ids,
    NumLib::LocalToGlobalIndexMap const& dof_table, std::size_t const mesh_id,
    int const variable_id, int const component_id,
    const std::vector<std::unique_ptr<ProcessLib::ParameterBase>>& parameters)
{
    DBUG("Constructing DirichletBoundaryCondition from config.");
    //! \ogs_file_param{boundary_condition__type}
    config.checkConfigParameter("type", "Dirichlet");

    //! \ogs_file_param{boundary_condition__Dirichlet__parameter}
    auto const param_name = config.getConfigParameter<std::string>("parameter");
    DBUG("Using parameter %s", param_name.c_str());

    auto& param = findParameter<double>(param_name, parameters, 1);

    return std::unique_ptr<DirichletBoundaryCondition>(
        new DirichletBoundaryCondition(param, std::move(mesh_node_ids),
                                       dof_table, mesh_id, variable_id,
                                       component_id));
}

}  // namespace ProcessLib
