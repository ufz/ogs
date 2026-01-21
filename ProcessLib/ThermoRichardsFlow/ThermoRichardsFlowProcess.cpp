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
    std::string name, MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    unsigned const integration_order,
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
        process_variables,
    ThermoRichardsFlowProcessData&& process_data,
    SecondaryVariableCollection&& secondary_variables,
    bool const use_monolithic_scheme, bool const is_linear)
    : Process(std::move(name), mesh, std::move(jacobian_assembler), parameters,
              integration_order, std::move(process_variables),
              std::move(secondary_variables), use_monolithic_scheme),
      AssemblyMixin<ThermoRichardsFlowProcess>{*_jacobian_assembler, is_linear,
                                               use_monolithic_scheme},
      _process_data(std::move(process_data))
{
    // For numerical Jacobian
    this->_jacobian_assembler->setNonDeformationComponentIDs({0, 1} /*T, p */);

    // TODO (naumov) remove ip suffix. Probably needs modification of the mesh
    // properties, s.t. there is no "overlapping" with cell/point data.
    // See getOrCreateMeshProperty.
    _integration_point_writer.emplace_back(
        std::make_unique<MeshLib::IntegrationPointWriter>(
            "saturation_ip", 1 /*n components*/, integration_order,
            local_assemblers_, &LocalAssemblerIF::getSaturation));

    _integration_point_writer.emplace_back(
        std::make_unique<MeshLib::IntegrationPointWriter>(
            "porosity_ip", 1 /*n components*/, integration_order,
            local_assemblers_, &LocalAssemblerIF::getPorosity));
}

void ThermoRichardsFlowProcess::initializeConcreteProcess(
    NumLib::LocalToGlobalIndexMap const& dof_table,
    MeshLib::Mesh const& mesh,
    unsigned const integration_order)
{
    ProcessLib::createLocalAssemblers<ThermoRichardsFlowLocalAssembler>(
        mesh.getDimension(), mesh.getElements(), dof_table, local_assemblers_,
        NumLib::IntegrationOrder{integration_order}, mesh.isAxiallySymmetric(),
        _process_data);

    auto add_secondary_variable = [&](std::string const& name,
                                      int const num_components,
                                      auto get_ip_values_function)
    {
        _secondary_variables.addSecondaryVariable(
            name,
            makeExtrapolator(num_components, getExtrapolator(),
                             local_assemblers_,
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
                               local_assemblers_);

    // Initialize local assemblers after all variables have been set.
    GlobalExecutor::executeMemberOnDereferenced(&LocalAssemblerIF::initialize,
                                                local_assemblers_,
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
        &LocalAssemblerIF::setInitialConditions, local_assemblers_,
        getDOFTables(x.size()), x, t, process_id);
}

std::vector<std::vector<std::string>>
ThermoRichardsFlowProcess::initializeAssemblyOnSubmeshes(
    std::vector<std::reference_wrapper<MeshLib::Mesh>> const& meshes)
{
    DBUG("TRF process initializeSubmeshOutput().");

    std::vector<std::vector<std::string>> per_process_residuum_names;

    if (_process_variables.size() == 1)  // monolithic
    {
        per_process_residuum_names = {{"HeatFlowRate", "MassFlowRate"}};
    }
    else  // staggered
    {
        per_process_residuum_names = {{"HeatFlowRate"}, {"MassFlowRate"}};
    }

    AssemblyMixin<ThermoRichardsFlowProcess>::initializeAssemblyOnSubmeshes(
        meshes, per_process_residuum_names);

    return per_process_residuum_names;
}

void ThermoRichardsFlowProcess::assembleConcreteProcess(
    const double t, double const dt, std::vector<GlobalVector*> const& x,
    std::vector<GlobalVector*> const& x_prev, int const process_id,
    GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b)
{
    DBUG("Assemble the equations for ThermoRichardsFlowProcess.");

    AssemblyMixin<ThermoRichardsFlowProcess>::assemble(t, dt, x, x_prev,
                                                       process_id, M, K, b);
}

void ThermoRichardsFlowProcess::assembleWithJacobianConcreteProcess(
    const double t, double const dt, std::vector<GlobalVector*> const& x,
    std::vector<GlobalVector*> const& x_prev, int const process_id,
    GlobalVector& b, GlobalMatrix& Jac)
{
    DBUG(
        "Assemble the Jacobian of ThermoRichardsFlow for the monolithic "
        "scheme.");

    AssemblyMixin<ThermoRichardsFlowProcess>::assembleWithJacobian(
        t, dt, x, x_prev, process_id, b, Jac);
}

void ThermoRichardsFlowProcess::preTimestepConcreteProcess(
    std::vector<GlobalVector*> const& /*x*/,
    const double /*t*/,
    const double /*dt*/,
    const int process_id)
{
    if (process_id == 0)
    {
        AssemblyMixin<ThermoRichardsFlowProcess>::updateActiveElements();
    }
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
        &LocalAssemblerIF::postTimestep, local_assemblers_,
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
        &LocalAssemblerIF::computeSecondaryVariable, local_assemblers_,
        getActiveElementIDs(), getDOFTables(x.size()), t, dt, x, x_prev,
        process_id);
}

}  // namespace ThermoRichardsFlow
}  // namespace ProcessLib
