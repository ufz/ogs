/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ThermoMechanicsProcess.h"

#include <cassert>

#include "MeshLib/Utils/IntegrationPointWriter.h"
#include "MeshLib/Utils/getOrCreateMeshProperty.h"
#include "NumLib/DOF/ComputeSparsityPattern.h"
#include "NumLib/DOF/DOFTableUtil.h"
#include "ProcessLib/Deformation/SolidMaterialInternalToSecondaryVariables.h"
#include "ProcessLib/SmallDeformation/CreateLocalAssemblers.h"
#include "ProcessLib/Utils/SetIPDataInitialConditions.h"
#include "ThermoMechanicsFEM.h"

namespace ProcessLib
{
namespace ThermoMechanics
{
template <int DisplacementDim>
ThermoMechanicsProcess<DisplacementDim>::ThermoMechanicsProcess(
    std::string name, MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    unsigned const integration_order,
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
        process_variables,
    ThermoMechanicsProcessData<DisplacementDim>&& process_data,
    SecondaryVariableCollection&& secondary_variables,
    bool const use_monolithic_scheme)
    : Process(std::move(name), mesh, std::move(jacobian_assembler), parameters,
              integration_order, std::move(process_variables),
              std::move(secondary_variables), use_monolithic_scheme),
      _process_data(std::move(process_data))
{
    _nodal_forces = MeshLib::getOrCreateMeshProperty<double>(
        mesh, "NodalForces", MeshLib::MeshItemType::Node, DisplacementDim);

    _heat_flux = MeshLib::getOrCreateMeshProperty<double>(
        mesh, "HeatFlowRate", MeshLib::MeshItemType::Node, 1);

    _integration_point_writer.emplace_back(
        std::make_unique<MeshLib::IntegrationPointWriter>(
            "sigma_ip",
            static_cast<int>(mesh.getDimension() == 2 ? 4 : 6) /*n components*/,
            integration_order, _local_assemblers,
            &LocalAssemblerInterface::getSigma));

    _integration_point_writer.emplace_back(
        std::make_unique<MeshLib::IntegrationPointWriter>(
            "epsilon_ip",
            static_cast<int>(mesh.getDimension() == 2 ? 4 : 6) /*n components*/,
            integration_order, _local_assemblers,
            &LocalAssemblerInterface::getEpsilon));

    _integration_point_writer.emplace_back(
        std::make_unique<MeshLib::IntegrationPointWriter>(
            "epsilon_m_ip",
            static_cast<int>(mesh.getDimension() == 2 ? 4 : 6) /*n components*/,
            integration_order, _local_assemblers,
            &LocalAssemblerInterface::getEpsilonMechanical));
}

template <int DisplacementDim>
bool ThermoMechanicsProcess<DisplacementDim>::isLinear() const
{
    return false;
}

template <int DisplacementDim>
MathLib::MatrixSpecifications
ThermoMechanicsProcess<DisplacementDim>::getMatrixSpecifications(
    const int process_id) const
{
    // For the monolithic scheme or the M process (deformation) in the staggered
    // scheme.
    if (_use_monolithic_scheme ||
        process_id == _process_data.mechanics_process_id)
    {
        auto const& l = *_local_to_global_index_map;
        return {l.dofSizeWithoutGhosts(), l.dofSizeWithoutGhosts(),
                &l.getGhostIndices(), &this->_sparsity_pattern};
    }

    // For staggered scheme and T process.
    auto const& l = *_local_to_global_index_map_single_component;
    return {l.dofSizeWithoutGhosts(), l.dofSizeWithoutGhosts(),
            &l.getGhostIndices(), &_sparsity_pattern_with_single_component};
}

// TODO [WW]: remove if (_use_monolithic_scheme) during the refactoring of the
// coupling part.
template <int DisplacementDim>
void ThermoMechanicsProcess<DisplacementDim>::constructDofTable()
{
    // Note: the heat conduction process and the mechanical process use the same
    // order of shape functions.

    if (_use_monolithic_scheme)
    {
        constructMonolithicProcessDofTable();
        return;
    }
    constructDofTableOfSpecifiedProcessStaggeredScheme(
        _process_data.mechanics_process_id);

    // TODO move the two data members somewhere else.
    // for extrapolation of secondary variables of stress or strain
    std::vector<MeshLib::MeshSubset> all_mesh_subsets_single_component{
        *_mesh_subset_all_nodes};
    _local_to_global_index_map_single_component.reset(
        new NumLib::LocalToGlobalIndexMap(
            std::move(all_mesh_subsets_single_component),
            // by location order is needed for output
            NumLib::ComponentOrder::BY_LOCATION));

    if (!_use_monolithic_scheme)
    {
        _sparsity_pattern_with_single_component =
            NumLib::computeSparsityPattern(
                *_local_to_global_index_map_single_component, _mesh);
    }
}

template <int DisplacementDim>
void ThermoMechanicsProcess<DisplacementDim>::initializeConcreteProcess(
    NumLib::LocalToGlobalIndexMap const& dof_table,
    MeshLib::Mesh const& mesh,
    unsigned const integration_order)
{
    ProcessLib::SmallDeformation::createLocalAssemblers<
        DisplacementDim, ThermoMechanicsLocalAssembler>(
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

    add_secondary_variable("sigma",
                           MathLib::KelvinVector::KelvinVectorType<
                               DisplacementDim>::RowsAtCompileTime,
                           &LocalAssemblerInterface::getIntPtSigma);

    add_secondary_variable("epsilon",
                           MathLib::KelvinVector::KelvinVectorType<
                               DisplacementDim>::RowsAtCompileTime,
                           &LocalAssemblerInterface::getIntPtEpsilon);

    //
    // enable output of internal variables defined by material models
    //
    ProcessLib::Deformation::solidMaterialInternalToSecondaryVariables<
        LocalAssemblerInterface>(_process_data.solid_materials,
                                 add_secondary_variable);

    setIPDataInitialConditions(_integration_point_writer, mesh.getProperties(),
                               _local_assemblers);

    // Initialize local assemblers after all variables have been set.
    GlobalExecutor::executeMemberOnDereferenced(
        &LocalAssemblerInterface::initialize, _local_assemblers,
        *_local_to_global_index_map);
}

template <int DisplacementDim>
void ThermoMechanicsProcess<DisplacementDim>::initializeBoundaryConditions(
    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const& media)
{
    if (_use_monolithic_scheme)
    {
        const int process_id_of_thermomechanics = 0;
        initializeProcessBoundaryConditionsAndSourceTerms(
            *_local_to_global_index_map, process_id_of_thermomechanics, media);
        return;
    }

    // Staggered scheme:
    // for the equations of heat conduction
    initializeProcessBoundaryConditionsAndSourceTerms(
        *_local_to_global_index_map_single_component,
        _process_data.heat_conduction_process_id, media);

    // for the equations of deformation.
    initializeProcessBoundaryConditionsAndSourceTerms(
        *_local_to_global_index_map, _process_data.mechanics_process_id, media);
}

template <int DisplacementDim>
void ThermoMechanicsProcess<DisplacementDim>::assembleConcreteProcess(
    const double t, double const dt, std::vector<GlobalVector*> const& x,
    std::vector<GlobalVector*> const& x_prev, int const process_id,
    GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b)
{
    DBUG("Assemble ThermoMechanicsProcess.");

    std::vector<NumLib::LocalToGlobalIndexMap const*> dof_table = {
        _local_to_global_index_map.get()};
    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    // Call global assembler for each local assembly item.
    GlobalExecutor::executeSelectedMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assemble, _local_assemblers,
        pv.getActiveElementIDs(), dof_table, t, dt, x, x_prev, process_id, M, K,
        b);
}

template <int DisplacementDim>
void ThermoMechanicsProcess<DisplacementDim>::
    assembleWithJacobianConcreteProcess(
        const double t, double const dt, std::vector<GlobalVector*> const& x,
        std::vector<GlobalVector*> const& x_prev, int const process_id,
        GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b, GlobalMatrix& Jac)
{
    DBUG("AssembleJacobian ThermoMechanicsProcess.");

    std::vector<NumLib::LocalToGlobalIndexMap const*> dof_tables;
    // For the monolithic scheme
    if (_use_monolithic_scheme)
    {
        DBUG(
            "Assemble the Jacobian of ThermoMechanics for the monolithic"
            " scheme.");
        dof_tables.emplace_back(_local_to_global_index_map.get());
    }
    else
    {
        // For the staggered scheme
        if (process_id == _process_data.heat_conduction_process_id)
        {
            DBUG(
                "Assemble the Jacobian equations of heat conduction process in "
                "ThermoMechanics for the staggered scheme.");
        }
        else
        {
            DBUG(
                "Assemble the Jacobian equations of mechanical process in "
                "ThermoMechanics for the staggered scheme.");
        }

        // For the flexible appearance order of processes in the coupling.
        if (_process_data.heat_conduction_process_id ==
            0)  // First: the heat conduction process
        {
            dof_tables.emplace_back(
                _local_to_global_index_map_single_component.get());
            dof_tables.emplace_back(_local_to_global_index_map.get());
        }
        else  // vice versa
        {
            dof_tables.emplace_back(_local_to_global_index_map.get());
            dof_tables.emplace_back(
                _local_to_global_index_map_single_component.get());
        }
    }

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    GlobalExecutor::executeSelectedMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assembleWithJacobian,
        _local_assemblers, pv.getActiveElementIDs(), dof_tables, t, dt, x,
        x_prev, process_id, M, K, b, Jac);

    // TODO (naumov): Refactor the copy rhs part. This is copy from HM.
    auto copyRhs = [&](int const variable_id, auto& output_vector)
    {
        if (_use_monolithic_scheme)
        {
            transformVariableFromGlobalVector(b, variable_id, *dof_tables[0],
                                              output_vector,
                                              std::negate<double>());
        }
        else
        {
            transformVariableFromGlobalVector(b, 0, *dof_tables[process_id],
                                              output_vector,
                                              std::negate<double>());
        }
    };
    if (_use_monolithic_scheme ||
        process_id == _process_data.heat_conduction_process_id)
    {
        copyRhs(0, *_heat_flux);
    }
    if (_use_monolithic_scheme ||
        process_id == _process_data.mechanics_process_id)
    {
        copyRhs(1, *_nodal_forces);
    }
}

template <int DisplacementDim>
void ThermoMechanicsProcess<DisplacementDim>::preTimestepConcreteProcess(
    std::vector<GlobalVector*> const& x, double const t, double const dt,
    const int process_id)
{
    DBUG("PreTimestep ThermoMechanicsProcess.");

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    assert(process_id < 2);

    if (process_id == _process_data.mechanics_process_id)
    {
        GlobalExecutor::executeSelectedMemberOnDereferenced(
            &LocalAssemblerInterface::preTimestep, _local_assemblers,
            pv.getActiveElementIDs(), *_local_to_global_index_map,
            *x[process_id], t, dt);
        return;
    }
}

template <int DisplacementDim>
void ThermoMechanicsProcess<DisplacementDim>::postTimestepConcreteProcess(
    std::vector<GlobalVector*> const& x,
    std::vector<GlobalVector*> const& x_prev, double const t, double const dt,
    int const process_id)
{
    if (process_id != 0)
    {
        return;
    }

    DBUG("PostTimestep ThermoMechanicsProcess.");
    std::vector<NumLib::LocalToGlobalIndexMap const*> dof_tables;
    auto const n_processes = x.size();
    dof_tables.reserve(n_processes);
    for (std::size_t process_id = 0; process_id < n_processes; ++process_id)
    {
        dof_tables.push_back(&getDOFTable(process_id));
    }

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];
    GlobalExecutor::executeSelectedMemberOnDereferenced(
        &LocalAssemblerInterface::postTimestep, _local_assemblers,
        pv.getActiveElementIDs(), dof_tables, x, x_prev, t, dt, process_id);
}

template <int DisplacementDim>
NumLib::LocalToGlobalIndexMap const&
ThermoMechanicsProcess<DisplacementDim>::getDOFTable(const int process_id) const
{
    if (_process_data.mechanics_process_id == process_id)
    {
        return *_local_to_global_index_map;
    }

    // For the equation of pressure
    return *_local_to_global_index_map_single_component;
}

template class ThermoMechanicsProcess<2>;
template class ThermoMechanicsProcess<3>;

}  // namespace ThermoMechanics
}  // namespace ProcessLib
