/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "GroundwaterFlowProcess.h"

#include <cassert>

#include "ProcessLib/Utils/CreateLocalAssemblers.h"

namespace ProcessLib
{
namespace GroundwaterFlow
{
GroundwaterFlowProcess::GroundwaterFlowProcess(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    std::vector<std::reference_wrapper<ProcessVariable>>&& process_variables,
    GroundwaterFlowProcessData&& process_data,
    SecondaryVariableCollection&& secondary_variables,
    NumLib::NamedFunctionCaller&& named_function_caller)
    : Process(mesh, std::move(jacobian_assembler), parameters,
              std::move(process_variables), std::move(secondary_variables),
              std::move(named_function_caller)),
      _process_data(std::move(process_data))
{
}

void GroundwaterFlowProcess::initializeConcreteProcess(
    NumLib::LocalToGlobalIndexMap const& dof_table,
    MeshLib::Mesh const& mesh,
    unsigned const integration_order)
{
    ProcessLib::createLocalAssemblers<LocalAssemblerData>(
        mesh.getDimension(), mesh.getElements(), dof_table, integration_order,
        _local_assemblers, _process_data);

    _secondary_variables.addSecondaryVariable(
        "darcy_velocity_x", 1,
        makeExtrapolator(
            getExtrapolator(), _local_assemblers,
            &GroundwaterFlowLocalAssemblerInterface::getIntPtDarcyVelocityX));

    if (mesh.getDimension() > 1) {
        _secondary_variables.addSecondaryVariable(
            "darcy_velocity_y", 1,
            makeExtrapolator(getExtrapolator(), _local_assemblers,
                             &GroundwaterFlowLocalAssemblerInterface::
                                 getIntPtDarcyVelocityY));
    }
    if (mesh.getDimension() > 2) {
        _secondary_variables.addSecondaryVariable(
            "darcy_velocity_z", 1,
            makeExtrapolator(getExtrapolator(), _local_assemblers,
                             &GroundwaterFlowLocalAssemblerInterface::
                                 getIntPtDarcyVelocityZ));
    }
}

void GroundwaterFlowProcess::assembleConcreteProcess(const double t,
                                                     GlobalVector const& x,
                                                     GlobalMatrix& M,
                                                     GlobalMatrix& K,
                                                     GlobalVector& b)
{
    DBUG("Assemble GroundwaterFlowProcess.");

    // Call global assembler for each local assembly item.
    GlobalExecutor::executeMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assemble, _local_assemblers,
        *_local_to_global_index_map, t, x, M, K, b);
}

void GroundwaterFlowProcess::assembleWithJacobianConcreteProcess(
    const double t, GlobalVector const& x, GlobalVector const& xdot,
    const double dxdot_dx, const double dx_dx, GlobalMatrix& M, GlobalMatrix& K,
    GlobalVector& b, GlobalMatrix& Jac)
{
    DBUG("AssembleWithJacobian GroundwaterFlowProcess.");

    // Call global assembler for each local assembly item.
    GlobalExecutor::executeMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assembleWithJacobian,
        _local_assemblers, *_local_to_global_index_map, t, x, xdot, dxdot_dx,
        dx_dx, M, K, b, Jac);
}

}   // namespace GroundwaterFlow
}   // namespace ProcessLib
