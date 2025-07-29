/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on 2025-07-29 11:42:36
 */

#include "TimeDecayDirichletBoundaryCondition.h"

#include "DirichletBoundaryConditionAuxiliaryFunctions.h"
#include "MathLib/LinAlg/LinAlg.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshSubset.h"
#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "NumLib/IndexValueVector.h"
#include "ParameterLib/Parameter.h"

namespace ProcessLib
{
std::unique_ptr<NumLib::LocalToGlobalIndexMap> createBoundaryDOFTable(
    int const variable_id, int const component_id,
    MeshLib::Mesh const& boundary_mesh,
    NumLib::LocalToGlobalIndexMap const& dof_table_bulk)
{
    checkParametersOfDirichletBoundaryCondition(boundary_mesh, dof_table_bulk,
                                                variable_id, component_id);

    std::vector<MeshLib::Node*> const& bc_nodes = boundary_mesh.getNodes();
    MeshLib::MeshSubset bc_mesh_subset(boundary_mesh, bc_nodes);

    // Create local DOF table from the BC mesh subset for the given variable
    // and component id.
    return dof_table_bulk.deriveBoundaryConstrainedMap(
        variable_id, {component_id}, std::move(bc_mesh_subset));
}

TimeDecayDirichletBoundaryCondition::TimeDecayDirichletBoundaryCondition(
    int const variable_id, int const component_id,
    MeshLib::Mesh const& boundary_mesh,
    NumLib::LocalToGlobalIndexMap const& dof_table_bulk,
    ParameterLib::Parameter<double> const& time_decay_parameter,
    double const lower_limit)
    : variable_id_(variable_id),
      component_id_(component_id),
      boundary_mesh_(boundary_mesh),
      dof_table_boundary_(createBoundaryDOFTable(
          variable_id, component_id, boundary_mesh, dof_table_bulk)),
      time_decay_parameter_(time_decay_parameter),
      lower_limit_(lower_limit)
{
}

void TimeDecayDirichletBoundaryCondition::setInitialValues(
    GlobalVector const& x) const
{
    initial_x_values_.clear();
    for (auto const& l : MeshLib::views::meshLocations(
             boundary_mesh_, MeshLib::MeshItemType::Node))
    {
        auto const global_index =
            dof_table_boundary_->getGlobalIndex(l, variable_id_, component_id_);
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
            initial_x_values_.push_back(x.get(global_index));
        }
    }
}

void TimeDecayDirichletBoundaryCondition::getEssentialBCValues(
    const double t, GlobalVector const& x,
    NumLib::IndexValueVector<GlobalIndexType>& bc_values) const
{
    DBUG("Apply TimeDecayDirichlet boundary bondition.");

    if (!initial_values_are_set_)
    {
        setInitialValues(x);
        initial_values_are_set_ = true;
        if (initial_x_values_.empty())
        {
            return;  // No node in the boundary.
        }
    }

    bc_values.ids.clear();
    bc_values.values.clear();

    bc_values.ids.reserve(initial_x_values_.size());
    bc_values.values.reserve(initial_x_values_.size());

    std::size_t idx = 0;
    for (auto const* node : boundary_mesh_.getNodes())
    {
        auto const node_id = node->getID();
        MeshLib::Location const l{boundary_mesh_.getID(),
                                  MeshLib::MeshItemType::Node, node_id};

        auto const global_index =
            dof_table_boundary_->getGlobalIndex(l, variable_id_, component_id_);

        if (global_index == NumLib::MeshComponentMap::nop)
        {
            continue;
        }
        // For the DDC approach (e.g. with PETSc option), the negative
        // index of global_index means that the entry by that index is a
        // ghost one, which should be dropped. Especially for PETSc
        // routines MatZeroRows and MatZeroRowsColumns, which are called
        // to apply the Dirichlet BC, the negative index is not accepted
        // like other matrix or vector PETSc routines. Therefore, the
        // following if-condition is applied.
        if (global_index >= 0) [[likely]]
        {
            ParameterLib::SpatialPosition pos;
            pos.setNodeID(node_id);
            pos.setCoordinates(*node);

            bc_values.values.push_back(
                time_decay_parameter_(t, pos).front() *
                    (initial_x_values_[idx] - lower_limit_) +
                lower_limit_);
            bc_values.ids.push_back(global_index);
            idx++;
        }
    }
}

}  // namespace ProcessLib
