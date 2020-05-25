/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
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
      source_term_mesh_id_(source_term_mesh_id),
      st_mesh_(st_mesh),
      variable_id_(variable_id),
      component_id_(component_id),
      parameter_(parameter)
{
    DBUG("Create NodalSourceTerm.");
}

void NodalSourceTerm::integrate(const double t, GlobalVector const& /*x*/,
                                GlobalVector& b, GlobalMatrix* /*jac*/) const
{
    DBUG("Assemble NodalSourceTerm.");

    for (MeshLib::Node const* const node : st_mesh_.getNodes())
    {
        auto const node_id = node->getID();
        MeshLib::Location const l{source_term_mesh_id_,
                                  MeshLib::MeshItemType::Node, node_id};
        auto const index = source_term_dof_table_->getGlobalIndex(
            l, variable_id_, component_id_);

        ParameterLib::SpatialPosition pos;
        pos.setNodeID(node_id);

        b.add(index, parameter_(t, pos).front());
    }
}

}  // namespace ProcessLib
