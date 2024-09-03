/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "WellboreSimulatorProcess.h"

#include <cassert>

#include "MeshLib/Utils/getOrCreateMeshProperty.h"
#include "NumLib/DOF/DOFTableUtil.h"
#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "ProcessLib/Utils/CreateLocalAssemblers.h"
#include "WellboreSimulatorFEM.h"

namespace ProcessLib
{
namespace WellboreSimulator
{
WellboreSimulatorProcess::WellboreSimulatorProcess(
    std::string name,
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    unsigned const integration_order,
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
        process_variables,
    WellboreSimulatorProcessData&& process_data,
    SecondaryVariableCollection&& secondary_variables)
    : Process(std::move(name), mesh, std::move(jacobian_assembler), parameters,
              integration_order, std::move(process_variables),
              std::move(secondary_variables)),
      _process_data(std::move(process_data))
{
}

void WellboreSimulatorProcess::initializeConcreteProcess(
    NumLib::LocalToGlobalIndexMap const& dof_table,
    MeshLib::Mesh const& mesh,
    unsigned const integration_order)
{
    ProcessLib::createLocalAssemblers<1, WellboreSimulatorFEM>(
        mesh.getElements(), dof_table, _local_assemblers,
        NumLib::IntegrationOrder{integration_order}, mesh.isAxiallySymmetric(),
        _process_data);

    auto add_secondary_variable = [&](std::string const& name,
                                      int const num_components,
                                      auto get_ip_values_function)
    {
        _secondary_variables.addSecondaryVariable(
            name,
            makeExtrapolator(num_components, getExtrapolator(),
                             _local_assemblers,
                             std::move(get_ip_values_function)));
    };

    add_secondary_variable(
        "vapor_mass_flow_rate", 1,
        &WellboreSimulatorLocalAssemblerInterface::getIntPtVaporMassFlowRate);

    add_secondary_variable(
        "liquid_mass_flow_rate", 1,
        &WellboreSimulatorLocalAssemblerInterface::getIntPtLiquidMassFlowRate);

    add_secondary_variable(
        "temperature", 1,
        &WellboreSimulatorLocalAssemblerInterface::getIntPtTemperature);

    add_secondary_variable(
        "dryness", 1,
        &WellboreSimulatorLocalAssemblerInterface::getIntPtDryness);

    add_secondary_variable(
        "vapor_volume_fraction", 1,
        &WellboreSimulatorLocalAssemblerInterface::getIntPtVaporVolumeFraction);

    _process_data.mesh_prop_density = MeshLib::getOrCreateMeshProperty<double>(
        const_cast<MeshLib::Mesh&>(mesh), "mix_density",
        MeshLib::MeshItemType::Cell, 1);
}

void WellboreSimulatorProcess::assembleConcreteProcess(
    const double t, double const dt, std::vector<GlobalVector*> const& x,
    std::vector<GlobalVector*> const& x_prev, int const process_id,
    GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b)
{
    DBUG("Assemble WellboreSimulator Process.");

    std::vector<NumLib::LocalToGlobalIndexMap const*> dof_tables;

    dof_tables.emplace_back(_local_to_global_index_map.get());

    // Call global assembler for each local assembly item.
    GlobalExecutor::executeSelectedMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assemble, _local_assemblers,
        getActiveElementIDs(), dof_tables, t, dt, x, x_prev, process_id, &M, &K,
        &b);
}

void WellboreSimulatorProcess::assembleWithJacobianConcreteProcess(
    const double t, double const dt, std::vector<GlobalVector*> const& x,
    std::vector<GlobalVector*> const& x_prev, int const process_id,
    GlobalVector& b, GlobalMatrix& Jac)
{
    DBUG("AssembleWithJacobian WellboreSimulator Process.");

    std::vector<NumLib::LocalToGlobalIndexMap const*> dof_tables;

    dof_tables.emplace_back(_local_to_global_index_map.get());

    // Call global assembler for each local assembly item.
    GlobalExecutor::executeSelectedMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assembleWithJacobian,
        _local_assemblers, getActiveElementIDs(), dof_tables, t, dt, x, x_prev,
        process_id, &b, &Jac);
}

void WellboreSimulatorProcess::computeSecondaryVariableConcrete(
    double const t,
    double const dt,
    std::vector<GlobalVector*> const& x,
    GlobalVector const& x_prev,
    int const process_id)
{
    std::vector<NumLib::LocalToGlobalIndexMap const*> dof_tables;
    dof_tables.reserve(x.size());
    dof_tables.push_back(_local_to_global_index_map.get());

    GlobalExecutor::executeSelectedMemberOnDereferenced(
        &WellboreSimulatorLocalAssemblerInterface::computeSecondaryVariable,
        _local_assemblers, getActiveElementIDs(), dof_tables, t, dt, x, x_prev,
        process_id);
}

void WellboreSimulatorProcess::postTimestepConcreteProcess(
    std::vector<GlobalVector*> const& x,
    std::vector<GlobalVector*> const& x_prev,
    const double t,
    const double dt,
    int const process_id)
{
    std::vector<NumLib::LocalToGlobalIndexMap const*> dof_tables;
    dof_tables.reserve(x.size());
    dof_tables.push_back(_local_to_global_index_map.get());

    GlobalExecutor::executeSelectedMemberOnDereferenced(
        &WellboreSimulatorLocalAssemblerInterface::postTimestep,
        _local_assemblers, getActiveElementIDs(), dof_tables, x, x_prev, t, dt,
        process_id);
}
}  // namespace WellboreSimulator
}  // namespace ProcessLib
