/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "SmallDeformationProcess.h"

#include <cassert>

#include "MeshLib/Utils/IntegrationPointWriter.h"
#include "MeshLib/Utils/getOrCreateMeshProperty.h"
#include "ProcessLib/Deformation/SolidMaterialInternalToSecondaryVariables.h"
#include "ProcessLib/Process.h"
#include "ProcessLib/Reflection/ReflectionForExtrapolation.h"
#include "ProcessLib/Reflection/ReflectionForIPWriters.h"
#include "ProcessLib/SmallDeformation/CreateLocalAssemblers.h"
#include "ProcessLib/Utils/SetIPDataInitialConditions.h"
#include "SmallDeformationFEM.h"

namespace ProcessLib
{
namespace SmallDeformation
{
template <int DisplacementDim>
SmallDeformationProcess<DisplacementDim>::SmallDeformationProcess(
    std::string name,
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    unsigned const integration_order,
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
        process_variables,
    SmallDeformationProcessData<DisplacementDim>&& process_data,
    SecondaryVariableCollection&& secondary_variables)
    : Process(std::move(name), mesh, std::move(jacobian_assembler), parameters,
              integration_order, std::move(process_variables),
              std::move(secondary_variables)),
      _process_data(std::move(process_data))
{
    _nodal_forces = MeshLib::getOrCreateMeshProperty<double>(
        mesh, "NodalForces", MeshLib::MeshItemType::Node, DisplacementDim);

    _material_forces = MeshLib::getOrCreateMeshProperty<double>(
        mesh, "MaterialForces", MeshLib::MeshItemType::Node, DisplacementDim);

    _process_data.principal_stress_vector[0] =
        MeshLib::getOrCreateMeshProperty<double>(
            const_cast<MeshLib::Mesh&>(mesh), "principal_stress_vector_1",
            MeshLib::MeshItemType::Cell, 3);

    _process_data.principal_stress_vector[1] =
        MeshLib::getOrCreateMeshProperty<double>(
            const_cast<MeshLib::Mesh&>(mesh), "principal_stress_vector_2",
            MeshLib::MeshItemType::Cell, 3);

    _process_data.principal_stress_vector[2] =
        MeshLib::getOrCreateMeshProperty<double>(
            const_cast<MeshLib::Mesh&>(mesh), "principal_stress_vector_3",
            MeshLib::MeshItemType::Cell, 3);

    _process_data.principal_stress_values =
        MeshLib::getOrCreateMeshProperty<double>(
            const_cast<MeshLib::Mesh&>(mesh), "principal_stress_values",
            MeshLib::MeshItemType::Cell, 3);

    ProcessLib::Reflection::addReflectedIntegrationPointWriters<
        DisplacementDim>(SmallDeformationLocalAssemblerInterface<
                             DisplacementDim>::getReflectionDataForOutput(),
                         _integration_point_writer, integration_order,
                         _local_assemblers);
}

template <int DisplacementDim>
bool SmallDeformationProcess<DisplacementDim>::isLinear() const
{
    return false;
}

template <int DisplacementDim>
void SmallDeformationProcess<DisplacementDim>::initializeConcreteProcess(
    NumLib::LocalToGlobalIndexMap const& dof_table,
    MeshLib::Mesh const& mesh,
    unsigned const integration_order)
{
    ProcessLib::SmallDeformation::createLocalAssemblers<
        DisplacementDim, SmallDeformationLocalAssembler>(
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

    ProcessLib::Reflection::addReflectedSecondaryVariables<DisplacementDim>(
        SmallDeformationLocalAssemblerInterface<
            DisplacementDim>::getReflectionDataForOutput(),
        _secondary_variables, getExtrapolator(), _local_assemblers);

    //
    // enable output of internal variables defined by material models
    //
    ProcessLib::Deformation::solidMaterialInternalToSecondaryVariables<
        LocalAssemblerInterface>(_process_data.solid_materials,
                                 add_secondary_variable);

    ProcessLib::Deformation::
        solidMaterialInternalVariablesToIntegrationPointWriter(
            _process_data.solid_materials, _local_assemblers,
            _integration_point_writer, integration_order);

    setIPDataInitialConditions(_integration_point_writer, mesh.getProperties(),
                               _local_assemblers);

    // Initialize local assemblers after all variables have been set.
    GlobalExecutor::executeMemberOnDereferenced(
        &LocalAssemblerInterface::initialize, _local_assemblers,
        *_local_to_global_index_map);
}

template <int DisplacementDim>
void SmallDeformationProcess<DisplacementDim>::assembleConcreteProcess(
    double const t, double const dt, std::vector<GlobalVector*> const& x,
    std::vector<GlobalVector*> const& x_prev, int const process_id,
    GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b)
{
    DBUG("Assemble SmallDeformationProcess.");

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_table = {std::ref(*_local_to_global_index_map)};
    // Call global assembler for each local assembly item.
    GlobalExecutor::executeSelectedMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assemble, _local_assemblers,
        pv.getActiveElementIDs(), dof_table, t, dt, x, x_prev, process_id, M, K,
        b);

    _global_output(t, process_id, M, K, b);
}

template <int DisplacementDim>
void SmallDeformationProcess<DisplacementDim>::
    assembleWithJacobianConcreteProcess(
        double const t, double const dt, std::vector<GlobalVector*> const& x,
        std::vector<GlobalVector*> const& x_prev, int const process_id,
        GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b, GlobalMatrix& Jac)
{
    DBUG("AssembleWithJacobian SmallDeformationProcess.");

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_table = {std::ref(*_local_to_global_index_map)};
    // Call global assembler for each local assembly item.
    GlobalExecutor::executeSelectedMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assembleWithJacobian,
        _local_assemblers, pv.getActiveElementIDs(), dof_table, t, dt, x,
        x_prev, process_id, M, K, b, Jac);

    transformVariableFromGlobalVector(b, 0, *_local_to_global_index_map,
                                      *_nodal_forces, std::negate<double>());

    _global_output(t, process_id, M, K, b, &Jac);
}

template <int DisplacementDim>
void SmallDeformationProcess<DisplacementDim>::postTimestepConcreteProcess(
    std::vector<GlobalVector*> const& x,
    std::vector<GlobalVector*> const& x_prev, double const t, double const dt,
    int const process_id)
{
    DBUG("PostTimestep SmallDeformationProcess.");
    std::vector<NumLib::LocalToGlobalIndexMap const*> dof_tables;
    dof_tables.reserve(x.size());
    std::generate_n(std::back_inserter(dof_tables), x.size(),
                    [&]() { return _local_to_global_index_map.get(); });

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];
    GlobalExecutor::executeSelectedMemberOnDereferenced(
        &LocalAssemblerInterface::postTimestep, _local_assemblers,
        pv.getActiveElementIDs(), dof_tables, x, x_prev, t, dt, process_id);

    std::unique_ptr<GlobalVector> material_forces;
    ProcessLib::SmallDeformation::writeMaterialForces(
        material_forces, _local_assemblers, *_local_to_global_index_map,
        *x[process_id]);

    material_forces->copyValues(*_material_forces);
}

template <int DisplacementDim>
void SmallDeformationProcess<DisplacementDim>::computeSecondaryVariableConcrete(
    double const t, double const dt, std::vector<GlobalVector*> const& x,
    GlobalVector const& x_prev, int const process_id)
{
    DBUG("Compute the secondary variables for SmallDeformationProcess.");
    std::vector<NumLib::LocalToGlobalIndexMap const*> dof_tables;
    dof_tables.reserve(x.size());
    std::generate_n(std::back_inserter(dof_tables), x.size(),
                    [&]() { return _local_to_global_index_map.get(); });

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];
    GlobalExecutor::executeSelectedMemberOnDereferenced(
        &LocalAssemblerInterface::computeSecondaryVariable, _local_assemblers,
        pv.getActiveElementIDs(), dof_tables, t, dt, x, x_prev, process_id);
}
template class SmallDeformationProcess<2>;
template class SmallDeformationProcess<3>;

}  // namespace SmallDeformation
}  // namespace ProcessLib
