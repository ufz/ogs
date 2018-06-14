/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "ProcessLib/BoundaryCondition/GenericNaturalBoundaryConditionLocalAssembler.h"
#include "ProcessLib/Utils/CreateLocalAssemblers.h"

namespace ProcessLib
{
namespace LIE
{
template <typename BoundaryConditionData,
          template <typename, typename, unsigned>
          class LocalAssemblerImplementation>
template <typename Data>
GenericNaturalBoundaryCondition<BoundaryConditionData,
                                LocalAssemblerImplementation>::
    GenericNaturalBoundaryCondition(
        typename std::enable_if<
            std::is_same<typename std::decay<BoundaryConditionData>::type,
                         typename std::decay<Data>::type>::value,
            bool>::type is_axially_symmetric,
        unsigned const integration_order, unsigned const shapefunction_order,
        NumLib::LocalToGlobalIndexMap const& dof_table_bulk,
        int const variable_id, int const component_id,
        unsigned const global_dim, MeshLib::Mesh const& bc_mesh, Data&& data,
        FractureProperty const& fracture_prop)
    : _data(std::forward<Data>(data)),
      _bc_mesh(bc_mesh),
      _integration_order(integration_order)
{
    assert(component_id < dof_table_bulk.getNumberOfComponents());

    std::vector<MeshLib::Node*> const& bc_nodes = _bc_mesh.getNodes();
    DBUG("Found %d nodes for Natural BCs for the variable %d and component %d",
         bc_nodes.size(), variable_id, component_id);

    MeshLib::MeshSubset bc_mesh_subset(_bc_mesh, bc_nodes);

    // Create local DOF table from the bc mesh subset for the given variable and
    // component id.
    _dof_table_boundary.reset(dof_table_bulk.deriveBoundaryConstrainedMap(
        variable_id, {component_id}, std::move(bc_mesh_subset)));

    createLocalAssemblers<LocalAssemblerImplementation>(
        global_dim, _bc_mesh.getElements(), *_dof_table_boundary,
        shapefunction_order, _local_assemblers, is_axially_symmetric,
        _integration_order, _data, fracture_prop, variable_id);
}

template <typename BoundaryConditionData,
          template <typename, typename, unsigned>
          class LocalAssemblerImplementation>
void GenericNaturalBoundaryCondition<
    BoundaryConditionData,
    LocalAssemblerImplementation>::applyNaturalBC(const double t,
                                                  const GlobalVector& x,
                                                  GlobalMatrix& K,
                                                  GlobalVector& b)
{
    GlobalExecutor::executeMemberOnDereferenced(
        &GenericNaturalBoundaryConditionLocalAssemblerInterface::assemble,
        _local_assemblers, *_dof_table_boundary, t, x, K, b);
}

}  // namespace LIE
}  // namespace ProcessLib
