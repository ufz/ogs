/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "HeatConductionProcess.h"

#include <cassert>

#include "MathLib/LinAlg/FinalizeMatrixAssembly.h"
#include "MathLib/LinAlg/FinalizeVectorAssembly.h"
#include "MeshLib/Utils/getOrCreateMeshProperty.h"
#include "NumLib/DOF/DOFTableUtil.h"
#include "ProcessLib/Utils/ComputeResiduum.h"
#include "ProcessLib/Utils/CreateLocalAssemblers.h"

namespace ProcessLib
{
namespace HeatConduction
{
HeatConductionProcess::HeatConductionProcess(
    std::string name,
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    unsigned const integration_order,
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
        process_variables,
    HeatConductionProcessData&& process_data,
    SecondaryVariableCollection&& secondary_variables)
    : Process(std::move(name), mesh, std::move(jacobian_assembler), parameters,
              integration_order, std::move(process_variables),
              std::move(secondary_variables)),
      _process_data(std::move(process_data))
{
    _heat_flux = MeshLib::getOrCreateMeshProperty<double>(
        mesh, "HeatFlowRate", MeshLib::MeshItemType::Node, 1);
}

void HeatConductionProcess::initializeConcreteProcess(
    NumLib::LocalToGlobalIndexMap const& dof_table,
    MeshLib::Mesh const& mesh,
    unsigned const integration_order)
{
    int const mesh_space_dimension = _process_data.mesh_space_dimension;

    ProcessLib::createLocalAssemblers<LocalAssemblerData>(
        mesh_space_dimension, mesh.getElements(), dof_table, _local_assemblers,
        NumLib::IntegrationOrder{integration_order}, mesh.isAxiallySymmetric(),
        _process_data);

    _secondary_variables.addSecondaryVariable(
        "heat_flux",
        makeExtrapolator(
            mesh_space_dimension, getExtrapolator(), _local_assemblers,
            &HeatConductionLocalAssemblerInterface::getIntPtHeatFlux));
}

void HeatConductionProcess::assembleConcreteProcess(
    const double t, double const dt, std::vector<GlobalVector*> const& x,
    std::vector<GlobalVector*> const& x_prev, int const process_id,
    GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b)
{
    DBUG("Assemble HeatConductionProcess.");

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_table = {std::ref(*_local_to_global_index_map)};
    // Call global assembler for each local assembly item.
    GlobalExecutor::executeSelectedMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assemble, _local_assemblers,
        pv.getActiveElementIDs(), dof_table, t, dt, x, x_prev, process_id, M, K,
        b);

    MathLib::finalizeMatrixAssembly(M);
    MathLib::finalizeMatrixAssembly(K);
    MathLib::finalizeVectorAssembly(b);

    auto const residuum = computeResiduum(dt, *x[0], *x_prev[0], M, K, b);
    transformVariableFromGlobalVector(residuum, 0 /*variable id*/,
                                      *_local_to_global_index_map, *_heat_flux,
                                      std::negate<double>());
}

void HeatConductionProcess::assembleWithJacobianConcreteProcess(
    const double t, double const dt, std::vector<GlobalVector*> const& x,
    std::vector<GlobalVector*> const& x_prev, int const process_id,
    GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b, GlobalMatrix& Jac)
{
    DBUG("AssembleWithJacobian HeatConductionProcess.");

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_table = {std::ref(*_local_to_global_index_map)};
    // Call global assembler for each local assembly item.
    GlobalExecutor::executeSelectedMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assembleWithJacobian,
        _local_assemblers, pv.getActiveElementIDs(), dof_table, t, dt, x,
        x_prev, process_id, M, K, b, Jac);

    transformVariableFromGlobalVector(b, 0 /*variable id*/,
                                      *_local_to_global_index_map, *_heat_flux,
                                      std::negate<double>());
}

void HeatConductionProcess::computeSecondaryVariableConcrete(
    double const t, double const dt, std::vector<GlobalVector*> const& x,
    GlobalVector const& x_prev, int const process_id)
{
    DBUG("Compute heat flux for HeatConductionProcess.");

    std::vector<NumLib::LocalToGlobalIndexMap const*> dof_tables;
    dof_tables.reserve(x.size());
    std::generate_n(std::back_inserter(dof_tables), x.size(),
                    [&]() { return _local_to_global_index_map.get(); });

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];
    GlobalExecutor::executeSelectedMemberOnDereferenced(
        &HeatConductionLocalAssemblerInterface::computeSecondaryVariable,
        _local_assemblers, pv.getActiveElementIDs(), dof_tables, t, dt, x,
        x_prev, process_id);
}

}  // namespace HeatConduction
}  // namespace ProcessLib
