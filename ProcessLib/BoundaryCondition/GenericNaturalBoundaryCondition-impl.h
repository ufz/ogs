/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "GenericNaturalBoundaryCondition.h"
#include "ProcessLib/Utils/CreateLocalAssemblers.h"
#include "MeshLib/MeshSearch/NodeSearch.h"
#include "GenericNaturalBoundaryConditionLocalAssembler.h"

namespace ProcessLib
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
            unsigned const>::type integration_order,
        NumLib::LocalToGlobalIndexMap const& dof_table_bulk,
        int const variable_id, int const component_id,
        unsigned const global_dim,
        std::vector<MeshLib::Element*>&& elements, Data&& data)
    : _data(std::forward<Data>(data)),
      _elements(std::move(elements)),
      _integration_order(integration_order)
{
    assert(component_id <
           static_cast<int>(dof_table_bulk.getNumberOfComponents()));

    std::vector<MeshLib::Node*> nodes = MeshLib::getUniqueNodes(_elements);

    auto const& mesh_subsets =
        dof_table_bulk.getMeshSubsets(variable_id, component_id);

    // TODO extend the node intersection to all parts of mesh_subsets, i.e.
    // to each of the MeshSubset in the mesh_subsets.
    _mesh_subset_all_nodes.reset(
        mesh_subsets.getMeshSubset(0).getIntersectionByNodes(nodes));
    std::unique_ptr<MeshLib::MeshSubsets> all_mesh_subsets{
        new MeshLib::MeshSubsets{_mesh_subset_all_nodes.get()}};

    // Create local DOF table from intersected mesh subsets for the given
    // variable and component ids.
    _dof_table_boundary.reset(dof_table_bulk.deriveBoundaryConstrainedMap(
        variable_id, component_id, std::move(all_mesh_subsets), _elements));

    createLocalAssemblers<LocalAssemblerImplementation>(
        global_dim, _elements, *_dof_table_boundary, _local_assemblers,
        false /* TODO */, _integration_order, _data);
}

template <typename BoundaryConditionData,
          template <typename, typename, unsigned>
          class LocalAssemblerImplementation>
GenericNaturalBoundaryCondition<
    BoundaryConditionData,
    LocalAssemblerImplementation>::~GenericNaturalBoundaryCondition()
{
    for (auto e : _elements)
        delete e;
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

}  // ProcessLib
