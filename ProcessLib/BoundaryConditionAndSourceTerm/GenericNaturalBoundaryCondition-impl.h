/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "GenericNaturalBoundaryConditionLocalAssembler.h"
#include "ProcessLib/BoundaryConditionAndSourceTerm/Utils/CreateLocalAssemblers.h"

namespace ProcessLib
{
template <typename BoundaryConditionData,
          template <typename /* shp fct */, int /* global dim */>
          class LocalAssemblerImplementation>
template <typename Data>
GenericNaturalBoundaryCondition<BoundaryConditionData,
                                LocalAssemblerImplementation>::
    GenericNaturalBoundaryCondition(
        unsigned const integration_order, unsigned const shapefunction_order,
        NumLib::LocalToGlobalIndexMap const& dof_table_bulk,
        int const variable_id, int const component_id,
        unsigned const global_dim, MeshLib::Mesh const& bc_mesh, Data&& data)
    : _data(std::forward<Data>(data)), _bc_mesh(bc_mesh)
{
    static_assert(std::is_same_v<typename std::decay_t<BoundaryConditionData>,
                                 typename std::decay_t<Data>>,
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

    if (!_bc_mesh.getProperties().template existsPropertyVector<std::size_t>(
            MeshLib::getBulkIDString(MeshLib::MeshItemType::Node)))
    {
        OGS_FATAL(
            "The required bulk node ids map does not exist in the boundary "
            "mesh '{:s}'.",
            _bc_mesh.getName());
    }

    std::vector<MeshLib::Node*> const& bc_nodes = _bc_mesh.getNodes();
    DBUG(
        "Found {:d} nodes for Natural BCs for the variable {:d} and component "
        "{:d}",
        bc_nodes.size(), variable_id, component_id);

    MeshLib::MeshSubset bc_mesh_subset(_bc_mesh, bc_nodes);

    // Create local DOF table from the BC mesh subset for the given variable and
    // component id.
    _dof_table_boundary = dof_table_bulk.deriveBoundaryConstrainedMap(
        variable_id, {component_id}, std::move(bc_mesh_subset));

    BoundaryConditionAndSourceTerm::createLocalAssemblers<
        LocalAssemblerImplementation>(
        global_dim, _bc_mesh.getElements(), *_dof_table_boundary,
        shapefunction_order, _local_assemblers,
        NumLib::IntegrationOrder{integration_order},
        _bc_mesh.isAxiallySymmetric(), _data);
}

template <typename BoundaryConditionData,
          template <typename /* shp fct */, int /* global dim */>
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
        _local_assemblers, *_dof_table_boundary, t, x, process_id, K, b, Jac);
}

}  // namespace ProcessLib
