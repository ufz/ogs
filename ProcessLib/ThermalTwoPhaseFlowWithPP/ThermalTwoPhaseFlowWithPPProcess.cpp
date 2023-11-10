/**
 * \file
 * \copyright
 * Copyright (c) 2012-2023, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ThermalTwoPhaseFlowWithPPProcess.h"

#include <cassert>

#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"
#include "MeshLib/PropertyVector.h"
#include "ProcessLib/Utils/CreateLocalAssemblers.h"
#include "ThermalTwoPhaseFlowWithPPLocalAssembler.h"

namespace ProcessLib
{
namespace ThermalTwoPhaseFlowWithPP
{
ThermalTwoPhaseFlowWithPPProcess::ThermalTwoPhaseFlowWithPPProcess(
    std::string name,
    MeshLib::Mesh& mesh,
    std::unique_ptr<AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    unsigned const integration_order,
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
        process_variables,
    ThermalTwoPhaseFlowWithPPProcessData&& process_data,
    SecondaryVariableCollection&& secondary_variables)
    : Process(std::move(name), mesh, std::move(jacobian_assembler), parameters,
              integration_order, std::move(process_variables),
              std::move(secondary_variables)),
      _process_data(std::move(process_data))
{
    DBUG("Create Nonisothermal TwoPhase Flow Process model.");
}

void ThermalTwoPhaseFlowWithPPProcess::initializeConcreteProcess(
    NumLib::LocalToGlobalIndexMap const& dof_table,
    MeshLib::Mesh const& mesh,
    unsigned const integration_order)
{
    ProcessLib::createLocalAssemblers<ThermalTwoPhaseFlowWithPPLocalAssembler>(
        mesh.getDimension(), mesh.getElements(), dof_table, _local_assemblers,
        NumLib::IntegrationOrder{integration_order}, mesh.isAxiallySymmetric(),
        _process_data);

    _secondary_variables.addSecondaryVariable(
        "saturation",
        makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                         &ThermalTwoPhaseFlowWithPPLocalAssemblerInterface::
                             getIntPtSaturation));

    _secondary_variables.addSecondaryVariable(
        "pressure_wetting",
        makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                         &ThermalTwoPhaseFlowWithPPLocalAssemblerInterface::
                             getIntPtWettingPressure));

    _secondary_variables.addSecondaryVariable(
        "liquid_molar_fraction_contaminant",
        makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                         &ThermalTwoPhaseFlowWithPPLocalAssemblerInterface::
                             getIntPtLiquidMolFracContaminant));

    _secondary_variables.addSecondaryVariable(
        "gas_molar_fraction_water",
        makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                         &ThermalTwoPhaseFlowWithPPLocalAssemblerInterface::
                             getIntPtGasMolFracWater));

    _secondary_variables.addSecondaryVariable(
        "gas_molar_fraction_contaminant",
        makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                         &ThermalTwoPhaseFlowWithPPLocalAssemblerInterface::
                             getIntPtGasMolFracContaminant));
}

void ThermalTwoPhaseFlowWithPPProcess::assembleConcreteProcess(
    const double t, double const dt, std::vector<GlobalVector*> const& x,
    std::vector<GlobalVector*> const& x_prev, int const process_id,
    GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b)
{
    DBUG("Assemble ThermalTwoPhaseFlowWithPPProcess.");

    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_table = {std::ref(*_local_to_global_index_map)};
    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    // Call global assembler for each local assembly item.
    GlobalExecutor::executeSelectedMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assemble, _local_assemblers,
        pv.getActiveElementIDs(), dof_table, t, dt, x, x_prev, process_id, M, K,
        b);
}

void ThermalTwoPhaseFlowWithPPProcess::assembleWithJacobianConcreteProcess(
    const double t, double const dt, std::vector<GlobalVector*> const& x,
    std::vector<GlobalVector*> const& x_prev, int const process_id,
    GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b, GlobalMatrix& Jac)
{
    DBUG("AssembleWithJacobian ThermalTwoPhaseFlowWithPPProcess.");

    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_table = {std::ref(*_local_to_global_index_map)};
    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    // Call global assembler for each local assembly item.
    GlobalExecutor::executeSelectedMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assembleWithJacobian,
        _local_assemblers, pv.getActiveElementIDs(), dof_table, t, dt, x,
        x_prev, process_id, M, K, b, Jac);
}
void ThermalTwoPhaseFlowWithPPProcess::preTimestepConcreteProcess(
    std::vector<GlobalVector*> const& x, double const t, double const delta_t,
    const int process_id)
{
    DBUG("PreTimestep ThermalTwoPhaseFlowWithPPProcess.");

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];
    GlobalExecutor::executeSelectedMemberOnDereferenced(
        &LocalAssemblerInterface::preTimestep, _local_assemblers,
        pv.getActiveElementIDs(), *_local_to_global_index_map, *x[process_id],
        t, delta_t);
}

}  // namespace ThermalTwoPhaseFlowWithPP
}  // namespace ProcessLib
