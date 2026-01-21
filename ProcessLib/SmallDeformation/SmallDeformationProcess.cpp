// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "SmallDeformationProcess.h"

#include <cassert>

#include "MeshLib/Utils/IntegrationPointWriter.h"
#include "MeshLib/Utils/getOrCreateMeshProperty.h"
#include "NumLib/Exceptions.h"
#include "ProcessLib/Deformation/SolidMaterialInternalToSecondaryVariables.h"
#include "ProcessLib/Output/CellAverageAlgorithm.h"
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
    std::string name, MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    unsigned const integration_order,
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
        process_variables,
    SmallDeformationProcessData<DisplacementDim>&& process_data,
    SecondaryVariableCollection&& secondary_variables, bool const is_linear)
    : Process(std::move(name), mesh, std::move(jacobian_assembler), parameters,
              integration_order, std::move(process_variables),
              std::move(secondary_variables)),
      AssemblyMixin<SmallDeformationProcess<DisplacementDim>>{
          *_jacobian_assembler, is_linear, true /* use_monolithic_scheme */},
      process_data_(std::move(process_data))
{
    // If numerical Jacobian assembler is used.
    if (this->_jacobian_assembler->isPerturbationEnabled())
    {
        OGS_FATAL(
            "Numerical Jacobian assembly is not supported for the "
            "SmallDeformationProcess.");
    }

    material_forces_ = MeshLib::getOrCreateMeshProperty<double>(
        mesh, "MaterialForces", MeshLib::MeshItemType::Node, DisplacementDim);

    process_data_.principal_stress_vector[0] =
        MeshLib::getOrCreateMeshProperty<double>(
            const_cast<MeshLib::Mesh&>(mesh), "principal_stress_vector_1",
            MeshLib::MeshItemType::Cell, 3);

    process_data_.principal_stress_vector[1] =
        MeshLib::getOrCreateMeshProperty<double>(
            const_cast<MeshLib::Mesh&>(mesh), "principal_stress_vector_2",
            MeshLib::MeshItemType::Cell, 3);

    process_data_.principal_stress_vector[2] =
        MeshLib::getOrCreateMeshProperty<double>(
            const_cast<MeshLib::Mesh&>(mesh), "principal_stress_vector_3",
            MeshLib::MeshItemType::Cell, 3);

    process_data_.principal_stress_values =
        MeshLib::getOrCreateMeshProperty<double>(
            const_cast<MeshLib::Mesh&>(mesh), "principal_stress_values",
            MeshLib::MeshItemType::Cell, 3);

    ProcessLib::Reflection::addReflectedIntegrationPointWriters<
        DisplacementDim>(SmallDeformationLocalAssemblerInterface<
                             DisplacementDim>::getReflectionDataForOutput(),
                         _integration_point_writer, integration_order,
                         local_assemblers_);
}

template <int DisplacementDim>
bool SmallDeformationProcess<DisplacementDim>::isLinear() const
{
    return AssemblyMixin<SmallDeformationProcess<DisplacementDim>>::isLinear();
}

template <int DisplacementDim>
void SmallDeformationProcess<DisplacementDim>::initializeConcreteProcess(
    NumLib::LocalToGlobalIndexMap const& dof_table,
    MeshLib::Mesh const& mesh,
    unsigned const integration_order)
{
    ProcessLib::SmallDeformation::createLocalAssemblers<
        DisplacementDim, SmallDeformationLocalAssembler>(
        mesh.getElements(), dof_table, local_assemblers_,
        NumLib::IntegrationOrder{integration_order}, mesh.isAxiallySymmetric(),
        process_data_);

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

    ProcessLib::Reflection::addReflectedSecondaryVariables<DisplacementDim>(
        SmallDeformationLocalAssemblerInterface<
            DisplacementDim>::getReflectionDataForOutput(),
        _secondary_variables, getExtrapolator(), local_assemblers_);

    //
    // enable output of internal variables defined by material models
    //
    ProcessLib::Deformation::solidMaterialInternalToSecondaryVariables<
        LocalAssemblerInterface>(process_data_.solid_materials,
                                 add_secondary_variable);

    ProcessLib::Deformation::
        solidMaterialInternalVariablesToIntegrationPointWriter(
            process_data_.solid_materials, local_assemblers_,
            _integration_point_writer, integration_order);

    setIPDataInitialConditions(_integration_point_writer, mesh.getProperties(),
                               local_assemblers_);

    // Initialize local assemblers after all variables have been set.
    GlobalExecutor::executeMemberOnDereferenced(
        &LocalAssemblerInterface::initialize, local_assemblers_,
        *_local_to_global_index_map);
}

template <int DisplacementDim>
void SmallDeformationProcess<DisplacementDim>::
    setInitialConditionsConcreteProcess(std::vector<GlobalVector*>& x,
                                        double const t, int const process_id)
{
    DBUG("Set initial conditions of SmallDeformationProcess.");

    GlobalExecutor::executeMemberOnDereferenced(
        &LocalAssemblerInterface::setInitialConditions, local_assemblers_,
        getDOFTables(x.size()), x, t, process_id);
}

template <int DisplacementDim>
void SmallDeformationProcess<DisplacementDim>::assembleConcreteProcess(
    double const t, double const dt, std::vector<GlobalVector*> const& x,
    std::vector<GlobalVector*> const& x_prev, int const process_id,
    GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b)
{
    DBUG("Assemble SmallDeformationProcess.");

    AssemblyMixin<SmallDeformationProcess<DisplacementDim>>::assemble(
        t, dt, x, x_prev, process_id, M, K, b);
}

template <int DisplacementDim>
void SmallDeformationProcess<DisplacementDim>::
    assembleWithJacobianConcreteProcess(
        double const t, double const dt, std::vector<GlobalVector*> const& x,
        std::vector<GlobalVector*> const& x_prev, int const process_id,
        GlobalVector& b, GlobalMatrix& Jac)
{
    DBUG("AssembleWithJacobian SmallDeformationProcess.");

    AssemblyMixin<SmallDeformationProcess<DisplacementDim>>::
        assembleWithJacobian(t, dt, x, x_prev, process_id, b, Jac);
}

template <int DisplacementDim>
void SmallDeformationProcess<DisplacementDim>::preTimestepConcreteProcess(
    std::vector<GlobalVector*> const& x, double const t, double const dt,
    const int process_id)
{
    DBUG("PreTimestep SmallDeformationProcess.");

    GlobalExecutor::executeSelectedMemberOnDereferenced(
        &LocalAssemblerInterface::preTimestep, local_assemblers_,
        getActiveElementIDs(), *_local_to_global_index_map, *x[process_id], t,
        dt);

    AssemblyMixin<
        SmallDeformationProcess<DisplacementDim>>::updateActiveElements();
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

    GlobalExecutor::executeSelectedMemberOnDereferenced(
        &LocalAssemblerInterface::postTimestep, local_assemblers_,
        getActiveElementIDs(), dof_tables, x, x_prev, t, dt, process_id);

    std::unique_ptr<GlobalVector> material_forces;
    ProcessLib::SmallDeformation::writeMaterialForces(
        material_forces, local_assemblers_, *_local_to_global_index_map,
        *x[process_id]);

    material_forces->copyValues(std::span{*material_forces_});
}

template <int DisplacementDim>
std::vector<std::vector<std::string>>
SmallDeformationProcess<DisplacementDim>::initializeAssemblyOnSubmeshes(
    std::vector<std::reference_wrapper<MeshLib::Mesh>> const& meshes)
{
    INFO("SmallDeformation process initializeSubmeshOutput().");
    std::vector<std::vector<std::string>> residuum_names{{"NodalForces"}};

    AssemblyMixin<SmallDeformationProcess<DisplacementDim>>::
        initializeAssemblyOnSubmeshes(meshes, residuum_names);

    return residuum_names;
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

    GlobalExecutor::executeSelectedMemberOnDereferenced(
        &LocalAssemblerInterface::computeSecondaryVariable, local_assemblers_,
        getActiveElementIDs(), dof_tables, t, dt, x, x_prev, process_id);

    computeCellAverages<DisplacementDim>(cell_average_data_, local_assemblers_);
}
template class SmallDeformationProcess<2>;
template class SmallDeformationProcess<3>;

}  // namespace SmallDeformation
}  // namespace ProcessLib
