/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ThermoRichardsFlowProcess.h"

#include <cassert>

#include "BaseLib/Error.h"
#include "MeshLib/Elements/Utils.h"
#include "MeshLib/Utils/getOrCreateMeshProperty.h"
#include "NumLib/DOF/ComputeSparsityPattern.h"
#include "NumLib/DOF/DOFTableUtil.h"
#include "ProcessLib/Process.h"
#include "ProcessLib/Utils/CreateLocalAssemblers.h"
#include "ProcessLib/Utils/SetIPDataInitialConditions.h"
#include "ThermoRichardsFlowFEM.h"

namespace ProcessLib
{
namespace ThermoRichardsFlow
{
ThermoRichardsFlowProcess::ThermoRichardsFlowProcess(
    std::string name,
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    unsigned const integration_order,
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
        process_variables,
    ThermoRichardsFlowProcessData&& process_data,
    SecondaryVariableCollection&& secondary_variables,
    bool const use_monolithic_scheme)
    : Process(std::move(name), mesh, std::move(jacobian_assembler), parameters,
              integration_order, std::move(process_variables),
              std::move(secondary_variables), use_monolithic_scheme),
      _process_data(std::move(process_data))
{
    _heat_flux = MeshLib::getOrCreateMeshProperty<double>(
        mesh, "HeatFlowRate", MeshLib::MeshItemType::Node, 1);

    _hydraulic_flow = MeshLib::getOrCreateMeshProperty<double>(
        mesh, "MassFlowRate", MeshLib::MeshItemType::Node, 1);

    // TODO (naumov) remove ip suffix. Probably needs modification of the mesh
    // properties, s.t. there is no "overlapping" with cell/point data.
    // See getOrCreateMeshProperty.
    _integration_point_writer.emplace_back(
        std::make_unique<MeshLib::IntegrationPointWriter>(
            "saturation_ip", 1 /*n components*/, integration_order,
            _local_assemblers, &LocalAssemblerIF::getSaturation));

    _integration_point_writer.emplace_back(
        std::make_unique<MeshLib::IntegrationPointWriter>(
            "porosity_ip", 1 /*n components*/, integration_order,
            _local_assemblers, &LocalAssemblerIF::getPorosity));
}

void ThermoRichardsFlowProcess::initializeConcreteProcess(
    NumLib::LocalToGlobalIndexMap const& dof_table,
    MeshLib::Mesh const& mesh,
    unsigned const integration_order)
{
    ProcessLib::createLocalAssemblers<ThermoRichardsFlowLocalAssembler>(
        mesh.getDimension(), mesh.getElements(), dof_table, _local_assemblers,
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

    add_secondary_variable("velocity", mesh.getDimension(),
                           &LocalAssemblerIF::getIntPtDarcyVelocity);

    add_secondary_variable("saturation", 1,
                           &LocalAssemblerIF::getIntPtSaturation);

    add_secondary_variable("porosity", 1, &LocalAssemblerIF::getIntPtPorosity);

    add_secondary_variable("dry_density_solid", 1,
                           &LocalAssemblerIF::getIntPtDryDensitySolid);

    _process_data.element_saturation = MeshLib::getOrCreateMeshProperty<double>(
        const_cast<MeshLib::Mesh&>(mesh), "saturation_avg",
        MeshLib::MeshItemType::Cell, 1);

    _process_data.element_porosity = MeshLib::getOrCreateMeshProperty<double>(
        const_cast<MeshLib::Mesh&>(mesh), "porosity_avg",
        MeshLib::MeshItemType::Cell, 1);

    setIPDataInitialConditions(_integration_point_writer, mesh.getProperties(),
                               _local_assemblers);

    // Initialize local assemblers after all variables have been set.
    GlobalExecutor::executeMemberOnDereferenced(&LocalAssemblerIF::initialize,
                                                _local_assemblers,
                                                *_local_to_global_index_map);
}

void ThermoRichardsFlowProcess::setInitialConditionsConcreteProcess(
    std::vector<GlobalVector*>& x, double const t, int const process_id)
{
    if (process_id != 0)
    {
        return;
    }
    DBUG("SetInitialConditions ThermoRichardsFlowProcess.");

    auto get_a_dof_table_func = [this](const int process_id) -> auto&
    {
        return getDOFTable(process_id);
    };
    GlobalExecutor::executeMemberOnDereferenced(
        &LocalAssemblerIF::setInitialConditions, _local_assemblers,
        NumLib::getDOFTables(x.size(), get_a_dof_table_func), x, t, process_id);
}

void ThermoRichardsFlowProcess::assembleConcreteProcess(
    const double t, double const dt, std::vector<GlobalVector*> const& x,
    std::vector<GlobalVector*> const& x_prev, int const process_id,
    GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b)
{
    DBUG("Assemble the equations for ThermoRichardsFlowProcess.");

    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_table = {std::ref(*_local_to_global_index_map)};
    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    // Call global assembler for each local assembly item.
    GlobalExecutor::executeSelectedMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assemble, _local_assemblers,
        pv.getActiveElementIDs(), dof_table, t, dt, x, x_prev, process_id, M, K,
        b);
}

void ThermoRichardsFlowProcess::assembleWithJacobianConcreteProcess(
    const double t, double const dt, std::vector<GlobalVector*> const& x,
    std::vector<GlobalVector*> const& x_prev, int const process_id,
    GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b, GlobalMatrix& Jac)
{
    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_tables;

    DBUG(
        "Assemble the Jacobian of ThermoRichardsFlow for the monolithic "
        "scheme.");
    dof_tables.emplace_back(*_local_to_global_index_map);

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    _pvma.assembleWithJacobian(_local_assemblers, pv.getActiveElementIDs(),
                               dof_tables, t, dt, x, x_prev, process_id, M, K,
                               b, Jac);

    auto copyRhs = [&](int const variable_id, auto& output_vector)
    {
        transformVariableFromGlobalVector(b, variable_id, dof_tables[0],
                                          output_vector, std::negate<double>());
    };

    copyRhs(0, *_heat_flux);
    copyRhs(1, *_hydraulic_flow);
}

void ThermoRichardsFlowProcess::postTimestepConcreteProcess(
    std::vector<GlobalVector*> const& x,
    std::vector<GlobalVector*> const& x_prev, double const t, double const dt,
    const int process_id)
{
    if (process_id != 0)
    {
        return;
    }

    DBUG("PostTimestep ThermoRichardsFlowProcess.");

    auto get_a_dof_table_func = [this](const int processe_id) -> auto&
    {
        return getDOFTable(processe_id);
    };
    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];
    GlobalExecutor::executeSelectedMemberOnDereferenced(
        &LocalAssemblerIF::postTimestep, _local_assemblers,
        pv.getActiveElementIDs(),
        NumLib::getDOFTables(x.size(), get_a_dof_table_func), x, x_prev, t, dt,
        process_id);
}

void ThermoRichardsFlowProcess::computeSecondaryVariableConcrete(
    const double t, const double dt, std::vector<GlobalVector*> const& x,
    GlobalVector const& x_prev, int const process_id)
{
    if (process_id != 0)
    {
        return;
    }
    DBUG(
        "Compute the secondary variables for "
        "ThermoRichardsFlowProcess.");

    auto get_a_dof_table_func = [this](const int processe_id) -> auto&
    {
        return getDOFTable(processe_id);
    };
    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    GlobalExecutor::executeSelectedMemberOnDereferenced(
        &LocalAssemblerIF::computeSecondaryVariable, _local_assemblers,
        pv.getActiveElementIDs(),
        NumLib::getDOFTables(x.size(), get_a_dof_table_func), t, dt, x, x_prev,
        process_id);
}

}  // namespace ThermoRichardsFlow
}  // namespace ProcessLib
