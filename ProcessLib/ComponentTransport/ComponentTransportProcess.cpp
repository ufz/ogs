/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ComponentTransportProcess.h"

#include <cassert>

#include "ProcessLib/Utils/CreateLocalAssemblers.h"

namespace ProcessLib
{
namespace ComponentTransport
{
ComponentTransportProcess::ComponentTransportProcess(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order,
    std::vector<std::reference_wrapper<ProcessVariable>>&& process_variables,
    ComponentTransportProcessData&& process_data,
    SecondaryVariableCollection&& secondary_variables,
    NumLib::NamedFunctionCaller&& named_function_caller)
    : Process(mesh, std::move(jacobian_assembler), parameters,
              integration_order, std::move(process_variables),
              std::move(secondary_variables), std::move(named_function_caller)),
      _process_data(std::move(process_data))
{
}

void ComponentTransportProcess::initializeConcreteProcess(
    NumLib::LocalToGlobalIndexMap const& dof_table,
    MeshLib::Mesh const& mesh,
    unsigned const integration_order)
{
    ProcessLib::ProcessVariable const& pv = getProcessVariables()[0];
    ProcessLib::createLocalAssemblers<LocalAssemblerData>(
        mesh.getDimension(), mesh.getElements(), dof_table,
        pv.getShapeFunctionOrder(), _local_assemblers,
        mesh.isAxiallySymmetric(), integration_order, _process_data);

    _secondary_variables.addSecondaryVariable(
        "darcy_velocity_x",
        makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                         &ComponentTransportLocalAssemblerInterface::
                             getIntPtDarcyVelocityX));

    if (mesh.getDimension() > 1)
    {
        _secondary_variables.addSecondaryVariable(
            "darcy_velocity_y",
            makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                             &ComponentTransportLocalAssemblerInterface::
                                 getIntPtDarcyVelocityY));
    }
    if (mesh.getDimension() > 2)
    {
        _secondary_variables.addSecondaryVariable(
            "darcy_velocity_z",
            makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                             &ComponentTransportLocalAssemblerInterface::
                                 getIntPtDarcyVelocityZ));
    }
}

void ComponentTransportProcess::assembleConcreteProcess(
    const double t,
    GlobalVector const& x,
    GlobalMatrix& M,
    GlobalMatrix& K,
    GlobalVector& b,
    StaggeredCouplingTerm const& coupling_term)
{
    DBUG("Assemble ComponentTransportProcess.");

    // Call global assembler for each local assembly item.
    GlobalExecutor::executeMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assemble, _local_assemblers,
        *_local_to_global_index_map, t, x, M, K, b, coupling_term);
}

void ComponentTransportProcess::assembleWithJacobianConcreteProcess(
    const double t, GlobalVector const& x, GlobalVector const& xdot,
    const double dxdot_dx, const double dx_dx, GlobalMatrix& M, GlobalMatrix& K,
    GlobalVector& b, GlobalMatrix& Jac,
    StaggeredCouplingTerm const& coupling_term)
{
    DBUG("AssembleWithJacobian ComponentTransportProcess.");

    // Call global assembler for each local assembly item.
    GlobalExecutor::executeMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assembleWithJacobian,
        _local_assemblers, *_local_to_global_index_map, t, x, xdot, dxdot_dx,
        dx_dx, M, K, b, Jac, coupling_term);
}

void ComponentTransportProcess::computeSecondaryVariableConcrete(
    double const t, GlobalVector const& x,
    StaggeredCouplingTerm const& coupling_term)
{
    DBUG("Compute the Darcy velocity for ComponentTransportProcess.");
    GlobalExecutor::executeMemberOnDereferenced(
        &ComponentTransportLocalAssemblerInterface::computeSecondaryVariable,
        _local_assemblers, *_local_to_global_index_map, t, x, coupling_term);
}

}  // namespace ComponentTransport
}  // namespace ProcessLib

