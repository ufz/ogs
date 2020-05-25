/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <numeric>

#include "MeshLib/MeshSearch/NodeSearch.h"
#include "ParameterLib/Utils.h"
#include "ProcessLib/Utils/CreateLocalAssemblers.h"

#include "NormalTractionBoundaryConditionLocalAssembler.h"

namespace ProcessLib
{
namespace NormalTractionBoundaryCondition
{
template <template <typename, typename, unsigned>
          class LocalAssemblerImplementation>
NormalTractionBoundaryCondition<LocalAssemblerImplementation>::
    NormalTractionBoundaryCondition(
        unsigned const integration_order, unsigned const shapefunction_order,
        NumLib::LocalToGlobalIndexMap const& dof_table_bulk,
        int const variable_id, unsigned const global_dim,
        MeshLib::Mesh const& bc_mesh,
        ParameterLib::Parameter<double> const& pressure)
    : bc_mesh_(bc_mesh),
      integration_order_(integration_order),
      pressure_(pressure)
{
    // Create component ids vector for the current variable.
    auto const& number_of_components =
        dof_table_bulk.getNumberOfVariableComponents(variable_id);
    std::vector<int> component_ids(number_of_components);
    std::iota(std::begin(component_ids), std::end(component_ids), 0);

    // BC mesh subset creation
    std::vector<MeshLib::Node*> const bc_nodes = bc_mesh_.getNodes();
    DBUG("Found {:d} nodes for Natural BCs for the variable {:d}",
         bc_nodes.size(), variable_id);

    MeshLib::MeshSubset bc_mesh_subset(bc_mesh_, bc_nodes);

    // Create local DOF table from the BC mesh subset for the given variable and
    // component ids.
    dof_table_boundary_.reset(dof_table_bulk.deriveBoundaryConstrainedMap(
        variable_id, component_ids, std::move(bc_mesh_subset)));

    createLocalAssemblers<LocalAssemblerImplementation>(
        global_dim, bc_mesh_.getElements(), *dof_table_boundary_,
        shapefunction_order, local_assemblers_, bc_mesh_.isAxiallySymmetric(),
        integration_order_, pressure_);
}

template <template <typename, typename, unsigned>
          class LocalAssemblerImplementation>
void NormalTractionBoundaryCondition<LocalAssemblerImplementation>::
    applyNaturalBC(const double t, std::vector<GlobalVector*> const& x,
                   int const /*process_id*/, GlobalMatrix& K, GlobalVector& b,
                   GlobalMatrix* Jac)
{
    GlobalExecutor::executeMemberOnDereferenced(
        &NormalTractionBoundaryConditionLocalAssemblerInterface::assemble,
        local_assemblers_, *dof_table_boundary_, t, x, K, b, Jac);
}

std::unique_ptr<NormalTractionBoundaryCondition<
    NormalTractionBoundaryConditionLocalAssembler>>
createNormalTractionBoundaryCondition(
    BaseLib::ConfigTree const& config, MeshLib::Mesh const& bc_mesh,
    NumLib::LocalToGlobalIndexMap const& dof_table, int const variable_id,
    unsigned const integration_order, unsigned const shapefunction_order,
    unsigned const global_dim,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters)
{
    DBUG("Constructing NormalTractionBoundaryCondition from config.");
    //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__type}
    config.checkConfigParameter("type", "NormalTraction");

    auto const parameter_name =
        //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__NormalTraction__parameter}
        config.getConfigParameter<std::string>("parameter");
    DBUG("Using parameter {:s}", parameter_name);

    auto const& pressure = ParameterLib::findParameter<double>(
        parameter_name, parameters, 1, &bc_mesh);
    return std::make_unique<NormalTractionBoundaryCondition<
        NormalTractionBoundaryConditionLocalAssembler>>(
        integration_order, shapefunction_order, dof_table, variable_id,
        global_dim, bc_mesh, pressure);
}

}  // namespace NormalTractionBoundaryCondition
}  // namespace ProcessLib
