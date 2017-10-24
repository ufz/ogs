/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "NodalSourceTerm.h"

#include <cassert>

namespace ProcessLib
{
NodalSourceTerm::NodalSourceTerm(const NumLib::LocalToGlobalIndexMap& dof_table,
                                 std::size_t const mesh_id,
                                 std::size_t const node_id,
                                 const int variable_id, const int component_id,
                                 double value)
    : _dof_table(dof_table),
      _mesh_id(mesh_id),
      _node_id(node_id),
      _variable_id(variable_id),
      _component_id(component_id),
      _value(value)
{
    DBUG("Create NodalSourceTerm.");
}

void NodalSourceTerm::integrateNodalSourceTerm(
    const double t,
    GlobalVector& b) const
{
    DBUG("Assemble NodalSourceTerm.");

    MeshLib::Location const l{_mesh_id, MeshLib::MeshItemType::Node, _node_id};
    auto const index =
        _dof_table.getGlobalIndex(l, _variable_id, _component_id);
    b.add(index, _value);
}

}  // namespace ProcessLib
