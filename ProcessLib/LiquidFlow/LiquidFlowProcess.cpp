/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file
 *
 * Created on August 19, 2016, 1:38 PM
 */

#include "LiquidFlowProcess.h"

#include <cassert>

#include "MeshLib/PropertyVector.h"

#include "ProcessLib/Utils/CreateLocalAssemblers.h"
#include "LiquidFlowLocalAssembler.h"
#include "LiquidFlowMaterialProperties.h"

#include "CreateLiquidFlowMaterialProperties.h"

namespace ProcessLib
{
namespace LiquidFlow
{
LiquidFlowProcess::LiquidFlowProcess(
    std::string name,
    MeshLib::Mesh& mesh,
    std::unique_ptr<AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    unsigned const integration_order,
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
        process_variables,
    SecondaryVariableCollection&& secondary_variables,
    MeshLib::PropertyVector<int> const* const material_ids,
    int const gravitational_axis_id,
    double const gravitational_acceleration,
    double const reference_temperature,
    BaseLib::ConfigTree const& config)
    : Process(std::move(name), mesh, std::move(jacobian_assembler), parameters,
              integration_order, std::move(process_variables),
              std::move(secondary_variables)),
      _gravitational_axis_id(gravitational_axis_id),
      _gravitational_acceleration(gravitational_acceleration),
      _reference_temperature(reference_temperature),
      _material_properties(
          createLiquidFlowMaterialProperties(config, parameters, material_ids))
{
    DBUG("Create Liquid flow process.");
}

void LiquidFlowProcess::initializeConcreteProcess(
    NumLib::LocalToGlobalIndexMap const& dof_table,
    MeshLib::Mesh const& mesh,
    unsigned const integration_order)
{
    const int process_id = 0;
    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];
    ProcessLib::createLocalAssemblers<LiquidFlowLocalAssembler>(
        mesh.getDimension(), mesh.getElements(), dof_table,
        pv.getShapeFunctionOrder(), _local_assemblers,
        mesh.isAxiallySymmetric(), integration_order, _gravitational_axis_id,
        _gravitational_acceleration, _reference_temperature,
        *_material_properties);

    _secondary_variables.addSecondaryVariable(
        "darcy_velocity",
        makeExtrapolator(
            mesh.getDimension(), getExtrapolator(), _local_assemblers,
            &LiquidFlowLocalAssemblerInterface::getIntPtDarcyVelocity));
}

void LiquidFlowProcess::assembleConcreteProcess(
    const double t, double const dt, std::vector<GlobalVector*> const& x,
    int const process_id, GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b)
{
    DBUG("Assemble LiquidFlowProcess.");

    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_table = {std::ref(*_local_to_global_index_map)};

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    // Call global assembler for each local assembly item.
    GlobalExecutor::executeSelectedMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assemble, _local_assemblers,
        pv.getActiveElementIDs(), dof_table, t, dt, x, process_id, M, K, b,
        _coupled_solutions);
}

void LiquidFlowProcess::assembleWithJacobianConcreteProcess(
    const double t, double const dt, std::vector<GlobalVector*> const& x,
    GlobalVector const& xdot, const double dxdot_dx, const double dx_dx,
    int const process_id, GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b,
    GlobalMatrix& Jac)
{
    DBUG("AssembleWithJacobian LiquidFlowProcess.");

    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_table = {std::ref(*_local_to_global_index_map)};
    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    // Call global assembler for each local assembly item.
    GlobalExecutor::executeSelectedMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assembleWithJacobian,
        _local_assemblers, pv.getActiveElementIDs(), dof_table, t, dt, x, xdot,
        dxdot_dx, dx_dx, process_id, M, K, b, Jac, _coupled_solutions);
}

void LiquidFlowProcess::computeSecondaryVariableConcrete(const double t,
                                                         GlobalVector const& x,
                                                         int const process_id)
{
    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    DBUG("Compute the velocity for LiquidFlowProcess.");
    GlobalExecutor::executeSelectedMemberOnDereferenced(
        &LiquidFlowLocalAssemblerInterface::computeSecondaryVariable,
        _local_assemblers, pv.getActiveElementIDs(),
        getDOFTable(process_id), t, x, _coupled_solutions);
}

}  // namespace LiquidFlow
}  // namespace ProcessLib
