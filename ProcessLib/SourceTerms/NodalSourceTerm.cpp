/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "NodalSourceTerm.h"

#include <cassert>

namespace ProcessLib
{
NodalSourceTerm::NodalSourceTerm(
    const NumLib::LocalToGlobalIndexMap& source_term_dof_table,
    std::size_t const bulk_mesh_id,
    MeshLib::Mesh const& st_mesh,
    const int variable_id,
    const int component_id,
    Parameter<double> const& parameter)
    : SourceTerm(source_term_dof_table),
      _bulk_mesh_id(bulk_mesh_id),
      _st_mesh(st_mesh),
      _variable_id(variable_id),
      _component_id(component_id),
      _parameter(parameter)
{
    DBUG("Create NodalSourceTerm.");
    if (!_st_mesh.getProperties().template existsPropertyVector<std::size_t>(
            "bulk_node_ids"))
    {
        OGS_FATAL(
            "Required mesh property \"bulk_node_ids\" does not exists on the "
            "source term mesh '%s'.", _st_mesh.getName().c_str());
    }
}

void NodalSourceTerm::integrate(const double t, GlobalVector& b) const
{
    DBUG("Assemble NodalSourceTerm.");

    auto const& bulk_node_ids_map =
        *_st_mesh.getProperties().template getPropertyVector<std::size_t>(
            "bulk_node_ids");
    for (MeshLib::Node const* const node : _st_mesh.getNodes())
    {
        auto const node_id = node->getID();
        MeshLib::Location const l{_bulk_mesh_id, MeshLib::MeshItemType::Node,
                                  bulk_node_ids_map[node_id]};
        auto const index = _source_term_dof_table.getGlobalIndex(
            l, _variable_id, _component_id);

        SpatialPosition pos;
        pos.setNodeID(node_id);

        b.add(index, _parameter(t, pos).front());
    }
}

}  // namespace ProcessLib
