/**
 * \file
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <numeric>

#include "MeshLib/MeshSearch/NodeSearch.h"
#include "ProcessLib/Utils/CreateLocalAssemblers.h"
#include "ProcessLib/Utils/ProcessUtils.h"

#include "NormalTractionBoundaryConditionLocalAssembler.h"

namespace ProcessLib
{
namespace NormalTractionBoundaryCondition
{
template <template <typename, typename, unsigned>
          class LocalAssemblerImplementation>
NormalTractionBoundaryCondition<LocalAssemblerImplementation>::
    NormalTractionBoundaryCondition(
        bool const is_axially_symmetric, unsigned const integration_order,
        unsigned const shapefunction_order,
        NumLib::LocalToGlobalIndexMap const& dof_table_bulk,
        int const variable_id, unsigned const global_dim,
        MeshLib::Mesh const& bc_mesh, Parameter<double> const& pressure)
    : _bc_mesh(bc_mesh),
      _integration_order(integration_order),
      _pressure(pressure)
{
    // Create component ids vector for the current variable.
    auto const& number_of_components =
        dof_table_bulk.getNumberOfVariableComponents(variable_id);
    std::vector<int> component_ids(number_of_components);
    std::iota(std::begin(component_ids), std::end(component_ids), 0);

    // BC mesh subset creation
    std::vector<MeshLib::Node*> const bc_nodes = _bc_mesh.getNodes();
    DBUG("Found %d nodes for Natural BCs for the variable %d", bc_nodes.size(),
         variable_id);

    MeshLib::MeshSubset bc_mesh_subset(_bc_mesh, bc_nodes);

    // Create local DOF table from the bc mesh subset for the given variable and
    // component ids.
    _dof_table_boundary.reset(dof_table_bulk.deriveBoundaryConstrainedMap(
        variable_id, component_ids, std::move(bc_mesh_subset)));

    createLocalAssemblers<LocalAssemblerImplementation>(
        global_dim, _bc_mesh.getElements(), *_dof_table_boundary,
        shapefunction_order, _local_assemblers, is_axially_symmetric,
        _integration_order, _pressure);
}

template <template <typename, typename, unsigned>
          class LocalAssemblerImplementation>
void NormalTractionBoundaryCondition<
    LocalAssemblerImplementation>::applyNaturalBC(const double t,
                                                  const GlobalVector& x,
                                                  GlobalMatrix& K,
                                                  GlobalVector& b)
{
    GlobalExecutor::executeMemberOnDereferenced(
        &NormalTractionBoundaryConditionLocalAssemblerInterface::assemble,
        _local_assemblers, *_dof_table_boundary, t, x, K, b);
}

std::unique_ptr<NormalTractionBoundaryCondition<
    NormalTractionBoundaryConditionLocalAssembler>>
createNormalTractionBoundaryCondition(
    BaseLib::ConfigTree const& config, MeshLib::Mesh const& bc_mesh,
    NumLib::LocalToGlobalIndexMap const& dof_table, int const variable_id,
    bool is_axially_symmetric, unsigned const integration_order,
    unsigned const shapefunction_order, unsigned const global_dim,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters)
{
    DBUG("Constructing NormalTractionBoundaryCondition from config.");
    //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__type}
    config.checkConfigParameter("type", "NormalTraction");

    auto const parameter_name =
        //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__NormalTraction__parameter}
        config.getConfigParameter<std::string>("parameter");
    DBUG("Using parameter %s", parameter_name.c_str());

    auto const& pressure = findParameter<double>(parameter_name, parameters, 1);
    return std::make_unique<NormalTractionBoundaryCondition<
        NormalTractionBoundaryConditionLocalAssembler>>(
        is_axially_symmetric, integration_order, shapefunction_order, dof_table,
        variable_id, global_dim, bc_mesh, pressure);
}

}  // namespace NormalTractionBoundaryCondition
}  // ProcessLib
