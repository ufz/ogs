/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on 2025-07-04 11:34:52
 */

#include "ReleaseNodalForce.h"

#include "MathLib/LinAlg/LinAlg.h"
#include "MeshLib/Mesh.h"
#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "ParameterLib/Parameter.h"

namespace ProcessLib
{
ReleaseNodalForce::ReleaseNodalForce(
    int const variable_id,
    MeshLib::Mesh const& boundary_mesh,
    std::unique_ptr<NumLib::LocalToGlobalIndexMap>& dof_table,
    ParameterLib::Parameter<double> const& time_decay_parameter)
    : variable_id_(variable_id),
      boundary_mesh_(boundary_mesh),
      dof_table_(std::move(dof_table)),
      time_decay_parameter_(time_decay_parameter)
{
}

void ReleaseNodalForce::set(GlobalVector const* r_neq)
{
    auto const& number_of_components =
        dof_table_->getNumberOfVariableComponents(variable_id_);

    if (!initial_release_nodal_force_.empty())
    {
        return;
    }

    initial_release_nodal_force_.reserve(boundary_mesh_.getNumberOfNodes() *
                                         number_of_components);

    // Iterate over all nodes in the boundary mesh and set the initial
    // released nodal.
    for (auto const* node : boundary_mesh_.getNodes())
    {
        auto const node_id = node->getID();
        MeshLib::Location const l{boundary_mesh_.getID(),
                                  MeshLib::MeshItemType::Node, node_id};
        for (int component_id = 0; component_id < number_of_components;
             ++component_id)
        {
            auto const global_index =
                dof_table_->getGlobalIndex(l, variable_id_, component_id);
            if (global_index == NumLib::MeshComponentMap::nop)
            {
                continue;
            }
            // For the DDC approach (e.g. with PETSc option), the negative
            // index of global_index means that the entry by that index is
            // a ghost one, which should be dropped. Especially for PETSc
            // routines MatZeroRows and MatZeroRowsColumns, which are
            // called to apply the Dirichlet BC, the negative index is not
            // accepted like other matrix or vector PETSc routines.
            // Therefore, the following if-condition is applied.
            if (global_index >= 0) [[likely]]
            {
                global_indices_.push_back(global_index);
                initial_release_nodal_force_.push_back(
                    r_neq->get(global_index));
                boundary_nodes_.push_back(node);
            }
        }
    }
}

void ReleaseNodalForce::applyNaturalBC(const double t,
                                       std::vector<GlobalVector*> const& /*x*/,
                                       int const /*process_id*/,
                                       GlobalMatrix* /*K*/, GlobalVector& b,
                                       GlobalMatrix* /*Jac*/)
{
    DBUG("Apply ReleaseNodalForce.");

    if (initial_release_nodal_force_.empty())
    {
        return;
    }

    std::vector<double> release_values;
    release_values.reserve(global_indices_.size());

    for (std::size_t i = 0; i < global_indices_.size(); ++i)
    {
        ParameterLib::SpatialPosition pos;
        auto const* node = boundary_nodes_[i];
        pos.setNodeID(node->getID());
        pos.setCoordinates(*node);

        release_values.push_back(
            -(1.0 - time_decay_parameter_(t, pos).front()) *
            initial_release_nodal_force_[i]);
    }

    b.add(global_indices_, release_values);
}

}  // namespace ProcessLib
