/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "RichardsComponentTransportProcess.h"

#include <cassert>

#include "ProcessLib/Utils/CreateLocalAssemblers.h"

namespace ProcessLib
{
namespace RichardsComponentTransport
{
RichardsComponentTransportProcess::RichardsComponentTransportProcess(
    std::string name,
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    unsigned const integration_order,
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
        process_variables,
    RichardsComponentTransportProcessData&& process_data,
    SecondaryVariableCollection&& secondary_variables,
    bool const use_monolithic_scheme)
    : Process(std::move(name), mesh, std::move(jacobian_assembler), parameters,
              integration_order, std::move(process_variables),
              std::move(secondary_variables), use_monolithic_scheme),
      process_data_(std::move(process_data))
{
}

void RichardsComponentTransportProcess::initializeConcreteProcess(
    NumLib::LocalToGlobalIndexMap const& dof_table,
    MeshLib::Mesh const& mesh,
    unsigned const integration_order)
{
    const int monolithic_process_id = 0;
    ProcessLib::ProcessVariable const& pv =
        getProcessVariables(monolithic_process_id)[0];
    ProcessLib::createLocalAssemblers<LocalAssemblerData>(
        mesh.getDimension(), mesh.getElements(), dof_table,
        pv.getShapeFunctionOrder(), local_assemblers_,
        mesh.isAxiallySymmetric(), integration_order, process_data_);

    secondary_variables_.addSecondaryVariable(
        "darcy_velocity",
        makeExtrapolator(mesh.getDimension(), getExtrapolator(),
                         local_assemblers_,
                         &RichardsComponentTransportLocalAssemblerInterface::
                             getIntPtDarcyVelocity));

    secondary_variables_.addSecondaryVariable(
        "saturation",
        makeExtrapolator(1, getExtrapolator(), local_assemblers_,
                         &RichardsComponentTransportLocalAssemblerInterface::
                             getIntPtSaturation));
}

void RichardsComponentTransportProcess::assembleConcreteProcess(
    const double t, double const dt, std::vector<GlobalVector*> const& x,
    std::vector<GlobalVector*> const& xdot, int const process_id,
    GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b)
{
    DBUG("Assemble RichardsComponentTransportProcess.");

    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
       dof_table = {std::ref(*local_to_global_index_map_)};
    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    // Call global assembler for each local assembly item.
    GlobalExecutor::executeSelectedMemberDereferenced(
        global_assembler_, &VectorMatrixAssembler::assemble, local_assemblers_,
        pv.getActiveElementIDs(), dof_table, t, dt, x, xdot, process_id, M, K,
        b, coupled_solutions_);
}

void RichardsComponentTransportProcess::assembleWithJacobianConcreteProcess(
    const double t, double const dt, std::vector<GlobalVector*> const& x,
    GlobalVector const& xdot, const double dxdot_dx, const double dx_dx,
    int const process_id, GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b,
    GlobalMatrix& Jac)
{
    DBUG("AssembleWithJacobian RichardsComponentTransportProcess.");

    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_table = {std::ref(*local_to_global_index_map_)};
    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    // Call global assembler for each local assembly item.
    GlobalExecutor::executeSelectedMemberDereferenced(
        global_assembler_, &VectorMatrixAssembler::assembleWithJacobian,
        local_assemblers_, pv.getActiveElementIDs(), dof_table, t, dt, x, xdot,
        dxdot_dx, dx_dx, process_id, M, K, b, Jac, coupled_solutions_);
}

}  // namespace RichardsComponentTransport
}  // namespace ProcessLib

