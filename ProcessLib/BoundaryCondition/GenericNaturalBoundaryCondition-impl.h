/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
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
            bool>::type is_axially_symmetric,
        unsigned const integration_order,
        unsigned const shapefunction_order,
        NumLib::LocalToGlobalIndexMap const& dof_table_bulk,
        int const variable_id, int const component_id,
        unsigned const global_dim, std::vector<MeshLib::Element*>&& elements,
        Data&& data)
    : _data(std::forward<Data>(data)),
      _elements(std::move(elements)),
      _integration_order(integration_order)
{
    // check basic data consistency
    if (variable_id >=
            static_cast<int>(dof_table_bulk.getNumberOfVariables()) ||
        component_id >=
            dof_table_bulk.getNumberOfVariableComponents(variable_id))
    {
        OGS_FATAL(
            "Variable id or component id too high. Actual values: (%d, %d), "
            "maximum values: (%d, %d).",
            variable_id, component_id, dof_table_bulk.getNumberOfVariables(),
            dof_table_bulk.getNumberOfVariableComponents(variable_id));
    }

    std::vector<MeshLib::Node*> nodes = MeshLib::getUniqueNodes(_elements);
    DBUG("Found %d nodes for Natural BCs for the variable %d and component %d",
         nodes.size(), variable_id, component_id);

    auto const& mesh_subset =
        dof_table_bulk.getMeshSubset(variable_id, component_id);

    MeshLib::MeshSubset bc_mesh_subset =
        mesh_subset.getIntersectionByNodes(nodes);

    // Create local DOF table from intersected mesh subsets for the given
    // variable and component ids.
    _dof_table_boundary.reset(dof_table_bulk.deriveBoundaryConstrainedMap(
        variable_id, {component_id}, std::move(bc_mesh_subset), _elements));

    createLocalAssemblers<LocalAssemblerImplementation>(
        global_dim, _elements, *_dof_table_boundary, shapefunction_order,
        _local_assemblers, is_axially_symmetric, _integration_order, _data);
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
