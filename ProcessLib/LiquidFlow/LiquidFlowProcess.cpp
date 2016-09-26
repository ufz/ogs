/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file   LiquidFlowProcess.cpp
 *
 * Created on August 19, 2016, 1:38 PM
 */

#include "LiquidFlowProcess.h"

#include <cassert>

#include "ProcessLib/Utils/CreateLocalAssemblers.h"
#include "LiquidFlowLocalAssembler.h"

namespace ProcessLib
{
namespace LiquidFlow
{
LiquidFlowProcess::LiquidFlowProcess(
    MeshLib::Mesh& mesh,
    std::unique_ptr<AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    std::vector<std::reference_wrapper<ProcessVariable>>&& process_variables,
    SecondaryVariableCollection&& secondary_variables,
    NumLib::NamedFunctionCaller&& named_function_caller,
    const bool compute_gravitational_term,
    BaseLib::ConfigTree const& config)
    : Process(mesh, std::move(jacobian_assembler), parameters,
              std::move(process_variables), std::move(secondary_variables),
              std::move(named_function_caller)),
      _compute_gravitational_term(compute_gravitational_term),
      _material_properties(LiquidFlowMaterialProperties(mesh, config))
{
    DBUG("Create Liquid flow process.");
}

void LiquidFlowProcess::initializeConcreteProcess(
    NumLib::LocalToGlobalIndexMap const& dof_table,
    MeshLib::Mesh const& mesh,
    unsigned const integration_order)
{
    ProcessLib::createLocalAssemblers<LiquidFlowLocalAssembler>(
        mesh.getDimension(), mesh.getElements(), dof_table, integration_order,
        _local_assemblers, _material_properties);

    _secondary_variables.addSecondaryVariable(
        "darcy_velocity_x", 1,
        makeExtrapolator(
            getExtrapolator(), _local_assemblers,
            &LiquidFlowLocalAssemblerInterface::getIntPtDarcyVelocityX));

    if (mesh.getDimension() > 1)
    {
        _secondary_variables.addSecondaryVariable(
            "darcy_velocity_y", 1,
            makeExtrapolator(
                getExtrapolator(), _local_assemblers,
                &LiquidFlowLocalAssemblerInterface::getIntPtDarcyVelocityY));
    }
    if (mesh.getDimension() > 2)
    {
        _secondary_variables.addSecondaryVariable(
            "darcy_velocity_z", 1,
            makeExtrapolator(
                getExtrapolator(), _local_assemblers,
                &LiquidFlowLocalAssemblerInterface::getIntPtDarcyVelocityZ));
    }
}

void LiquidFlowProcess::assembleConcreteProcess(const double t,
                                                GlobalVector const& p,
                                                GlobalMatrix& M,
                                                GlobalMatrix& K,
                                                GlobalVector& b)
{
    DBUG("Assemble LiquidFlowProcess.");
    // Call global assembler for each local assembly item.
    GlobalExecutor::executeMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assemble, _local_assemblers,
        *_local_to_global_index_map, t, p, M, K, b);
}

void LiquidFlowProcess::assembleWithJacobianConcreteProcess(
    const double t, GlobalVector const& x, GlobalVector const& xdot,
    const double dxdot_dx, const double dx_dx, GlobalMatrix& M, GlobalMatrix& K,
    GlobalVector& b, GlobalMatrix& Jac)
{
    DBUG("AssembleWithJacobian LiquidFlowProcess.");

    // Call global assembler for each local assembly item.
    GlobalExecutor::executeMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assembleWithJacobian,
        _local_assemblers, *_local_to_global_index_map, t, x, xdot, dxdot_dx,
        dx_dx, M, K, b, Jac);
}

void LiquidFlowProcess::preTimestepConcreteProcess(GlobalVector const& x,
                                                   const double /*t*/,
                                                   const double /*delta_t*/)
{
    DBUG("New time step");

    _p_previous_timestep =
        MathLib::MatrixVectorTraits<GlobalVector>::newInstance(x);
}

void LiquidFlowProcess::preIterationConcreteProcess(const unsigned /*iter*/,
                                                    GlobalVector const& /*x*/)
{
    // TODO
}
}  // end of namespace
}  // end of namespace
