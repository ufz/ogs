/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "GenericNonuniformNaturalBoundaryCondition.h"
#include "GenericNonuniformNaturalBoundaryConditionLocalAssembler.h"
#include "MeshLib/MeshSearch/NodeSearch.h"
#include "ProcessLib/Utils/CreateLocalAssemblers.h"

namespace ProcessLib
{
template <typename BoundaryConditionData,
          template <typename, typename, unsigned>
          class LocalAssemblerImplementation>
template <typename Data>
GenericNonuniformNaturalBoundaryCondition<BoundaryConditionData,
                                          LocalAssemblerImplementation>::
    GenericNonuniformNaturalBoundaryCondition(
        unsigned const integration_order, unsigned const shapefunction_order,
        unsigned const global_dim,
        std::unique_ptr<MeshLib::Mesh>&& boundary_mesh, Data&& data)
    : _data(std::forward<Data>(data)), _boundary_mesh(std::move(boundary_mesh))
{
    static_assert(std::is_same<typename std::decay<BoundaryConditionData>::type,
                               typename std::decay<Data>::type>::value,
                  "Type mismatch between declared and passed BC data.");

    // check basic data consistency
    if (_data.variable_id_bulk >=
            static_cast<int>(_data.dof_table_bulk.getNumberOfVariables()) ||
        _data.component_id_bulk >=
            _data.dof_table_bulk.getNumberOfVariableComponents(
                _data.variable_id_bulk))
    {
        OGS_FATAL(
            "Variable id or component id too high. Actual values: (%d, %d), "
            "maximum values: (%d, %d).",
            _data.variable_id_bulk, _data.component_id_bulk,
            _data.dof_table_bulk.getNumberOfVariables(),
            _data.dof_table_bulk.getNumberOfVariableComponents(
                _data.variable_id_bulk));
    }

    if (_boundary_mesh->getDimension() + 1 != global_dim)
    {
        OGS_FATAL(
            "The dimension of the given boundary mesh (%d) is not by one lower "
            "than the bulk dimension (%d).",
            _boundary_mesh->getDimension(), global_dim);
    }

    constructDofTable();

    createLocalAssemblers<LocalAssemblerImplementation>(
        global_dim, _boundary_mesh->getElements(), *_dof_table_boundary,
        shapefunction_order, _local_assemblers,
        _boundary_mesh->isAxiallySymmetric(), integration_order, _data);
}

template <typename BoundaryConditionData,
          template <typename, typename, unsigned>
          class LocalAssemblerImplementation>
void GenericNonuniformNaturalBoundaryCondition<
    BoundaryConditionData, LocalAssemblerImplementation>::constructDofTable()
{
    // construct one-component DOF-table for the surface mesh
    _mesh_subset_all_nodes.reset(
        new MeshLib::MeshSubset(*_boundary_mesh, _boundary_mesh->getNodes()));

    std::vector<MeshLib::MeshSubset> all_mesh_subsets{*_mesh_subset_all_nodes};

    std::vector<int> vec_var_n_components{1};

    _dof_table_boundary = std::make_unique<NumLib::LocalToGlobalIndexMap>(
        std::move(all_mesh_subsets), vec_var_n_components,
        NumLib::ComponentOrder::BY_LOCATION);
}

template <typename BoundaryConditionData,
          template <typename, typename, unsigned>
          class LocalAssemblerImplementation>
void GenericNonuniformNaturalBoundaryCondition<
    BoundaryConditionData,
    LocalAssemblerImplementation>::applyNaturalBC(const double t,
                                                  const GlobalVector& x,
                                                  GlobalMatrix& K,
                                                  GlobalVector& b,
                                                  GlobalMatrix* Jac)
{
    GlobalExecutor::executeMemberOnDereferenced(
        &GenericNonuniformNaturalBoundaryConditionLocalAssemblerInterface::
            assemble,
        _local_assemblers, *_dof_table_boundary, t, x, K, b, Jac);
}

}  // namespace ProcessLib
