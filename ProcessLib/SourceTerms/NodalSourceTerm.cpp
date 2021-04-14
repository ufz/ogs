/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
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
    std::unique_ptr<NumLib::LocalToGlobalIndexMap> source_term_dof_table,
    std::size_t const source_term_mesh_id,
    MeshLib::Mesh const& st_mesh,
    const int variable_id,
    const int component_id,
    ParameterLib::Parameter<double> const& parameter)
    : SourceTerm(std::move(source_term_dof_table)),
      _source_term_mesh_id(source_term_mesh_id),
      _st_mesh(st_mesh),
      _variable_id(variable_id),
      _component_id(component_id),
      _parameter(parameter)
{
    DBUG("Create NodalSourceTerm.");
}

void NodalSourceTerm::integrate(const double t, GlobalVector const& /*x*/,
                                GlobalVector& b, GlobalMatrix* /*jac*/) const
{
    DBUG("Assemble NodalSourceTerm.");

    for (MeshLib::Node const* const node : _st_mesh.getNodes())
    {
        auto const node_id = node->getID();
        MeshLib::Location const l{_source_term_mesh_id,
                                  MeshLib::MeshItemType::Node, node_id};
        auto const global_index = _source_term_dof_table->getGlobalIndex(
            l, _variable_id, _component_id);

        if (global_index == NumLib::MeshComponentMap::nop)
        {
            continue;
        }
        // For the DDC approach (e.g. with PETSc option), the negative
        // index of global_index means that the entry by that index is a ghost
        // one, which should be dropped. Especially for PETSc routines
        // MatZeroRows and MatZeroRowsColumns, which are called to apply the
        // Dirichlet BC, the negative index is not accepted like other matrix or
        // vector PETSc routines. Therefore, the following if-condition is
        // applied.
        if (global_index >= 0)
        {
            ParameterLib::SpatialPosition pos;
            pos.setNodeID(node_id);
            pos.setCoordinates(*node);

            b.add(global_index, _parameter(t, pos).front());
        }
    }
}

}  // namespace ProcessLib
