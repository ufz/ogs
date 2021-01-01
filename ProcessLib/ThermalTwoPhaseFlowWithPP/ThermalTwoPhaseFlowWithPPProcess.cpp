/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ThermalTwoPhaseFlowWithPPProcess.h"

#include <cassert>

#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"
#include "MeshLib/PropertyVector.h"
#include "MaterialLib/PorousMedium/Porosity/Porosity.h"
#include "MaterialLib/PorousMedium/Storage/Storage.h"
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
    SecondaryVariableCollection&& secondary_variables,
    BaseLib::ConfigTree const& /*config*/,
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
    /*curves*/)
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
    const int monolithic_process_id = 0;
    ProcessLib::ProcessVariable const& pv =
        getProcessVariables(monolithic_process_id)[0];
    ProcessLib::createLocalAssemblers<ThermalTwoPhaseFlowWithPPLocalAssembler>(
        mesh.getDimension(), mesh.getElements(), dof_table,
        pv.getShapeFunctionOrder(), _local_assemblers,
        mesh.isAxiallySymmetric(), integration_order, _process_data);

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
}

void ThermalTwoPhaseFlowWithPPProcess::assembleConcreteProcess(
    const double t, double const dt, std::vector<GlobalVector*> const& x,
    std::vector<GlobalVector*> const& xdot, int const process_id,
    GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b)
{
    DBUG("Assemble ThermalTwoPhaseFlowWithPPProcess.");

    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
       dof_table = {std::ref(*_local_to_global_index_map)};
    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    // Call global assembler for each local assembly item.
    GlobalExecutor::executeSelectedMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assemble, _local_assemblers,
        pv.getActiveElementIDs(), dof_table, t, dt, x, xdot, process_id, M, K,
        b);
}

void ThermalTwoPhaseFlowWithPPProcess::assembleWithJacobianConcreteProcess(
    const double t, double const dt, std::vector<GlobalVector*> const& x,
    std::vector<GlobalVector*> const& xdot, const double dxdot_dx,
    const double dx_dx, int const process_id, GlobalMatrix& M, GlobalMatrix& K,
    GlobalVector& b, GlobalMatrix& Jac)
{
    DBUG("AssembleWithJacobian ThermalTwoPhaseFlowWithPPProcess.");

    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
       dof_table = {std::ref(*_local_to_global_index_map)};
    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    // Call global assembler for each local assembly item.
    GlobalExecutor::executeSelectedMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assembleWithJacobian,
        _local_assemblers, pv.getActiveElementIDs(), dof_table, t, dt, x, xdot,
        dxdot_dx, dx_dx, process_id, M, K, b, Jac);
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
