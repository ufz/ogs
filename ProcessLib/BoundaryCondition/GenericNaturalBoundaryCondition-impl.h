/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "GenericNaturalBoundaryConditionLocalAssembler.h"
#include "ProcessLib/Utils/CreateLocalAssemblers.h"

namespace ProcessLib
{
template <typename BoundaryConditionData,
          template <typename, typename, unsigned>
          class LocalAssemblerImplementation>
template <typename Data>
GenericNaturalBoundaryCondition<BoundaryConditionData,
                                LocalAssemblerImplementation>::
    GenericNaturalBoundaryCondition(
        unsigned const integration_order, unsigned const shapefunction_order,
        NumLib::LocalToGlobalIndexMap const& dof_table_bulk,
        int const variable_id, int const component_id,
        unsigned const global_dim, MeshLib::Mesh const& bc_mesh, Data&& data)
    : data_(std::forward<Data>(data)), bc_mesh_(bc_mesh)
{
    static_assert(std::is_same<typename std::decay<BoundaryConditionData>::type,
                               typename std::decay<Data>::type>::value,
                  "Type mismatch between declared and passed BC data.");

    // check basic data consistency
    if (variable_id >=
            static_cast<int>(dof_table_bulk.getNumberOfVariables()) ||
        component_id >=
            dof_table_bulk.getNumberOfVariableComponents(variable_id))
    {
        OGS_FATAL(
            "Variable id or component id too high. Actual values: ({:d}, "
            "{:d}), "
            "maximum values: ({:d}, {:d}).",
            variable_id, component_id, dof_table_bulk.getNumberOfVariables(),
            dof_table_bulk.getNumberOfVariableComponents(variable_id));
    }

    if (!bc_mesh_.getProperties().template existsPropertyVector<std::size_t>(
            "bulk_node_ids"))
    {
        OGS_FATAL(
            "The required bulk node ids map does not exist in the boundary "
            "mesh '{:s}'.",
            bc_mesh_.getName());
    }

    std::vector<MeshLib::Node*> const& bc_nodes = bc_mesh_.getNodes();
    DBUG(
        "Found {:d} nodes for Natural BCs for the variable {:d} and component "
        "{:d}",
        bc_nodes.size(), variable_id, component_id);

    MeshLib::MeshSubset bc_mesh_subset(bc_mesh_, bc_nodes);

    // Create local DOF table from the BC mesh subset for the given variable and
    // component id.
    dof_table_boundary_.reset(dof_table_bulk.deriveBoundaryConstrainedMap(
        variable_id, {component_id}, std::move(bc_mesh_subset)));

    createLocalAssemblers<LocalAssemblerImplementation>(
        global_dim, bc_mesh_.getElements(), *dof_table_boundary_,
        shapefunction_order, local_assemblers_, bc_mesh_.isAxiallySymmetric(),
        integration_order, data_);
}

template <typename BoundaryConditionData,
          template <typename, typename, unsigned>
          class LocalAssemblerImplementation>
void GenericNaturalBoundaryCondition<BoundaryConditionData,
                                     LocalAssemblerImplementation>::
    applyNaturalBC(const double t,
                   std::vector<GlobalVector*> const& x,
                   int const process_id,
                   GlobalMatrix& K,
                   GlobalVector& b,
                   GlobalMatrix* Jac)
{
    GlobalExecutor::executeMemberOnDereferenced(
        &GenericNaturalBoundaryConditionLocalAssemblerInterface::assemble,
        local_assemblers_, *dof_table_boundary_, t, x, process_id, K, b, Jac);
}

}  // namespace ProcessLib
