// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "ThermoRichardsFlowProcess.h"

#include <cassert>

#include "BaseLib/Error.h"
#include "MeshLib/Elements/Utils.h"
#include "MeshLib/Utils/getOrCreateMeshProperty.h"
#include "NumLib/DOF/ComputeSparsityPattern.h"
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
    // For numerical Jacobian
    this->_jacobian_assembler->setNonDeformationComponentIDs({0, 1} /*T, p */);

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

    GlobalExecutor::executeMemberOnDereferenced(
        &LocalAssemblerIF::setInitialConditions, _local_assemblers,
        getDOFTables(x.size()), x, t, process_id);
}

void ThermoRichardsFlowProcess::assembleConcreteProcess(
    const double t, double const dt, std::vector<GlobalVector*> const& x,
    std::vector<GlobalVector*> const& x_prev, int const process_id,
    GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b)
{
    DBUG("Assemble the equations for ThermoRichardsFlowProcess.");

    std::vector<NumLib::LocalToGlobalIndexMap const*> dof_table = {
        _local_to_global_index_map.get()};

    // Call global assembler for each local assembly item.
    GlobalExecutor::executeSelectedMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assemble, _local_assemblers,
        getActiveElementIDs(), dof_table, t, dt, x, x_prev, process_id, &M, &K,
        &b);
}

void ThermoRichardsFlowProcess::assembleWithJacobianConcreteProcess(
    const double t, double const dt, std::vector<GlobalVector*> const& x,
    std::vector<GlobalVector*> const& x_prev, int const process_id,
    GlobalVector& b, GlobalMatrix& Jac)
{
    std::vector<NumLib::LocalToGlobalIndexMap const*> dof_tables;

    DBUG(
        "Assemble the Jacobian of ThermoRichardsFlow for the monolithic "
        "scheme.");
    dof_tables.emplace_back(_local_to_global_index_map.get());

    _pvma.assembleWithJacobian(_local_assemblers, getActiveElementIDs(),
                               dof_tables, t, dt, x, x_prev, process_id, b,
                               Jac);

    auto copyRhs = [&](int const variable_id, auto& output_vector)
    {
        transformVariableFromGlobalVector(b, variable_id, *dof_tables[0],
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

    GlobalExecutor::executeSelectedMemberOnDereferenced(
        &LocalAssemblerIF::postTimestep, _local_assemblers,
        getActiveElementIDs(), getDOFTables(x.size()), x, x_prev, t, dt,
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
    DBUG("Compute the secondary variables for ThermoRichardsFlowProcess.");

    GlobalExecutor::executeSelectedMemberOnDereferenced(
        &LocalAssemblerIF::computeSecondaryVariable, _local_assemblers,
        getActiveElementIDs(), getDOFTables(x.size()), t, dt, x, x_prev,
        process_id);
}

}  // namespace ThermoRichardsFlow
}  // namespace ProcessLib
