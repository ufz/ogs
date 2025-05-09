/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "TwoPhaseFlowWithPPProcess.h"

#include "ProcessLib/Utils/CreateLocalAssemblers.h"
#include "TwoPhaseFlowWithPPLocalAssembler.h"

namespace ProcessLib
{
namespace TwoPhaseFlowWithPP
{
TwoPhaseFlowWithPPProcess::TwoPhaseFlowWithPPProcess(
    std::string name,
    MeshLib::Mesh& mesh,
    std::unique_ptr<AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    unsigned const integration_order,
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
        process_variables,
    TwoPhaseFlowWithPPProcessData&& process_data,
    SecondaryVariableCollection&& secondary_variables)
    : Process(std::move(name), mesh, std::move(jacobian_assembler), parameters,
              integration_order, std::move(process_variables),
              std::move(secondary_variables)),
      _process_data(std::move(process_data))
{
    DBUG("Create TwoPhaseFlowProcess with PP model.");
}

void TwoPhaseFlowWithPPProcess::setInitialConditionsConcreteProcess(
    std::vector<GlobalVector*>& x, double const t, int const process_id)
{
    DBUG("SetInitialConditions ThermoRichardsMechanicsProcess.");

    GlobalExecutor::executeMemberOnDereferenced(
        &TwoPhaseFlowWithPPLocalAssemblerInterface::setInitialConditions,
        _local_assemblers, getDOFTables(x.size()), x, t, process_id);
}

void TwoPhaseFlowWithPPProcess::initializeConcreteProcess(
    NumLib::LocalToGlobalIndexMap const& dof_table,
    MeshLib::Mesh const& mesh,
    unsigned const integration_order)
{
    ProcessLib::createLocalAssemblers<TwoPhaseFlowWithPPLocalAssembler>(
        mesh.getDimension(), mesh.getElements(), dof_table, _local_assemblers,
        NumLib::IntegrationOrder{integration_order}, mesh.isAxiallySymmetric(),
        _process_data);

    _secondary_variables.addSecondaryVariable(
        "saturation",
        makeExtrapolator(
            1, getExtrapolator(), _local_assemblers,
            &TwoPhaseFlowWithPPLocalAssemblerInterface::getIntPtSaturation));

    _secondary_variables.addSecondaryVariable(
        "pressure_wet",
        makeExtrapolator(
            1, getExtrapolator(), _local_assemblers,
            &TwoPhaseFlowWithPPLocalAssemblerInterface::getIntPtWetPressure));
}

void TwoPhaseFlowWithPPProcess::assembleConcreteProcess(
    const double t, double const dt, std::vector<GlobalVector*> const& x,
    std::vector<GlobalVector*> const& x_prev, int const process_id,
    GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b)
{
    DBUG("Assemble TwoPhaseFlowWithPPProcess.");

    std::vector<NumLib::LocalToGlobalIndexMap const*> dof_table = {
        _local_to_global_index_map.get()};

    // Call global assembler for each local assembly item.
    GlobalExecutor::executeSelectedMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assemble, _local_assemblers,
        getActiveElementIDs(), dof_table, t, dt, x, x_prev, process_id, &M, &K,
        &b);
}

void TwoPhaseFlowWithPPProcess::assembleWithJacobianConcreteProcess(
    const double t, double const dt, std::vector<GlobalVector*> const& x,
    std::vector<GlobalVector*> const& x_prev, int const process_id,
    GlobalVector& b, GlobalMatrix& Jac)
{
    DBUG("AssembleWithJacobian TwoPhaseFlowWithPPProcess.");

    std::vector<NumLib::LocalToGlobalIndexMap const*> dof_table = {
        _local_to_global_index_map.get()};

    // Call global assembler for each local assembly item.
    GlobalExecutor::executeSelectedMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assembleWithJacobian,
        _local_assemblers, getActiveElementIDs(), dof_table, t, dt, x, x_prev,
        process_id, &b, &Jac);
}

}  // namespace TwoPhaseFlowWithPP
}  // namespace ProcessLib
