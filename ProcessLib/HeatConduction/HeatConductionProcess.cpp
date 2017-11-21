/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "HeatConductionProcess.h"

#include <cassert>

#include "ProcessLib/Utils/CreateLocalAssemblers.h"

namespace ProcessLib
{
namespace HeatConduction
{
HeatConductionProcess::HeatConductionProcess(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order,
    std::vector<std::reference_wrapper<ProcessVariable>>&& process_variables,
    HeatConductionProcessData&& process_data,
    SecondaryVariableCollection&& secondary_variables,
    NumLib::NamedFunctionCaller&& named_function_caller)
    : Process(mesh, std::move(jacobian_assembler), parameters,
              integration_order, std::move(process_variables),
              std::move(secondary_variables), std::move(named_function_caller)),
      _process_data(std::move(process_data))
{
}

void HeatConductionProcess::preTimestepConcreteProcess(GlobalVector const& x,
                                            const double /*t*/,
                                            const double /*delta_t*/,
                                            const int /*process_id*/)
{
    if (!_x_previous_timestep)
    {
        _x_previous_timestep =
            MathLib::MatrixVectorTraits<GlobalVector>::newInstance(x);
    }
    else
    {
        auto& x0 = *_x_previous_timestep;
        MathLib::LinAlg::copy(x, x0);
    }
}

void HeatConductionProcess::initializeConcreteProcess(
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
        "heat_flux_x",
        makeExtrapolator(
            1, getExtrapolator(), _local_assemblers,
            &HeatConductionLocalAssemblerInterface::getIntPtHeatFluxX));

    if (mesh.getDimension() > 1)
    {
        _secondary_variables.addSecondaryVariable(
            "heat_flux_y",
            makeExtrapolator(
                1, getExtrapolator(), _local_assemblers,
                &HeatConductionLocalAssemblerInterface::getIntPtHeatFluxY));
    }
    if (mesh.getDimension() > 2)
    {
        _secondary_variables.addSecondaryVariable(
            "heat_flux_z",
            makeExtrapolator(
                1, getExtrapolator(), _local_assemblers,
                &HeatConductionLocalAssemblerInterface::getIntPtHeatFluxZ));
    }
}

void HeatConductionProcess::assembleConcreteProcess(const double t,
                                                    GlobalVector const& x,
                                                    GlobalMatrix& M,
                                                    GlobalMatrix& K,
                                                    GlobalVector& b)
{
    DBUG("Assemble HeatConductionProcess.");

    // Call global assembler for each local assembly item.
    GlobalExecutor::executeMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assemble, _local_assemblers,
        *_local_to_global_index_map, t, x, M, K, b, _coupled_solutions);
}

void HeatConductionProcess::assembleWithJacobianConcreteProcess(
    const double t, GlobalVector const& x, GlobalVector const& xdot,
    const double dxdot_dx, const double dx_dx, GlobalMatrix& M, GlobalMatrix& K,
    GlobalVector& b, GlobalMatrix& Jac)
{
    DBUG("AssembleWithJacobian HeatConductionProcess.");

    // Call global assembler for each local assembly item.
    GlobalExecutor::executeMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assembleWithJacobian,
        _local_assemblers, *_local_to_global_index_map, t, x, xdot, dxdot_dx,
        dx_dx, M, K, b, Jac, _coupled_solutions);
}

void HeatConductionProcess::computeSecondaryVariableConcrete(
    const double t, GlobalVector const& x)
{
    DBUG("Compute heat flux for HeatConductionProcess.");
    GlobalExecutor::executeMemberOnDereferenced(
            &HeatConductionLocalAssemblerInterface::computeSecondaryVariable,
            _local_assemblers, *_local_to_global_index_map, t, x,
            _coupled_solutions);
}

}  // namespace HeatConduction
}  // namespace ProcessLib
