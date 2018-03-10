/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <cassert>

#include "BaseLib/Functional.h"
#include "MeshLib/Elements/Utils.h"
#include "NumLib/DOF/ComputeSparsityPattern.h"
#include "ProcessLib/SmallDeformation/CreateLocalAssemblers.h"
#include "ProcessLib/Process.h"

#include "ThermoMechanicalPhaseFieldFEM.h"
#include "ThermoMechanicalPhaseFieldProcessData.h"
#include "ThermoMechanicalPhaseFieldProcess.h"

namespace ProcessLib
{
namespace ThermoMechanicalPhaseField
{
template <int DisplacementDim>
ThermoMechanicalPhaseFieldProcess<DisplacementDim>::ThermoMechanicalPhaseFieldProcess(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order,
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
        process_variables,
    ThermoMechanicalPhaseFieldProcessData<DisplacementDim>&& process_data,
    SecondaryVariableCollection&& secondary_variables,
    NumLib::NamedFunctionCaller&& named_function_caller,
    bool const use_monolithic_scheme,
    int const mechanics_related_process_id,
    int const phase_field_process_id,
    int const heat_conduction_process_id)
    : Process(mesh, std::move(jacobian_assembler), parameters,
              integration_order, std::move(process_variables),
              std::move(secondary_variables), std::move(named_function_caller),
              use_monolithic_scheme),
      _process_data(std::move(process_data)),
      _mechanics_related_process_id(mechanics_related_process_id),
      _phase_field_process_id(phase_field_process_id),
      _heat_conduction_process_id(heat_conduction_process_id)
{
}

template <int DisplacementDim>
bool ThermoMechanicalPhaseFieldProcess<DisplacementDim>::isLinear() const
{
    return false;
}

template <int DisplacementDim>
MathLib::MatrixSpecifications
ThermoMechanicalPhaseFieldProcess<DisplacementDim>::getMatrixSpecifications(
    const int process_id) const
{
    if (process_id == _mechanics_related_process_id)
    {
        auto const& l = *_local_to_global_index_map;
        return {l.dofSizeWithoutGhosts(), l.dofSizeWithoutGhosts(),
                &l.getGhostIndices(), &this->_sparsity_pattern};
    }

    // For staggered scheme and phase field process or heat conduction.
    auto const& l = *_local_to_global_index_map_single_component;
    return {l.dofSizeWithoutGhosts(), l.dofSizeWithoutGhosts(),
            &l.getGhostIndices(), &_sparsity_pattern_with_single_component};
}

template <int DisplacementDim>
NumLib::LocalToGlobalIndexMap const&
ThermoMechanicalPhaseFieldProcess<DisplacementDim>::getDOFTable(const int process_id) const
{
    if (process_id == _mechanics_related_process_id)
    {
        return *_local_to_global_index_map;
    }

    // For the equation of phasefield or heat conduction.
    return *_local_to_global_index_map_single_component;
}

template <int DisplacementDim>
NumLib::LocalToGlobalIndexMap&
ThermoMechanicalPhaseFieldProcess<DisplacementDim>::getDOFTableByProcessID(
    const int process_id) const
{
    if (process_id == _mechanics_related_process_id)
    {
        return *_local_to_global_index_map;
    }

    // For the equation of phasefield or heat conduction.
    return *_local_to_global_index_map_single_component;
}

template <int DisplacementDim>
void ThermoMechanicalPhaseFieldProcess<DisplacementDim>::constructDofTable()
{
    // Create single component dof in every of the mesh's nodes.
    _mesh_subset_all_nodes =
        std::make_unique<MeshLib::MeshSubset>(_mesh, &_mesh.getNodes());

    // TODO move the two data members somewhere else.
    // for extrapolation of secondary variables of stress or strain
    std::vector<MeshLib::MeshSubsets> all_mesh_subsets_single_component;
    all_mesh_subsets_single_component.emplace_back(
        _mesh_subset_all_nodes.get());
    _local_to_global_index_map_single_component =
        std::make_unique<NumLib::LocalToGlobalIndexMap>(
            std::move(all_mesh_subsets_single_component),
            // by location order is needed for output
            NumLib::ComponentOrder::BY_LOCATION);

    assert(_local_to_global_index_map_single_component);

    if (_use_monolithic_scheme)
    {
        std::vector<MeshLib::MeshSubsets> all_mesh_subsets;
        all_mesh_subsets.emplace_back(_mesh_subset_all_nodes.get());
        all_mesh_subsets.emplace_back(_mesh_subset_all_nodes.get());

        const int monolithic_process_id =
            0;  // Only one process in the monolithic scheme.
        std::generate_n(
            std::back_inserter(all_mesh_subsets),
            getProcessVariables(monolithic_process_id)[2]
                .get()
                .getNumberOfComponents(),
            [&]() {
                return MeshLib::MeshSubsets{_mesh_subset_all_nodes.get()};
            });

        std::vector<int> const vec_n_components{1, 1, DisplacementDim};
        _local_to_global_index_map =
            std::make_unique<NumLib::LocalToGlobalIndexMap>(
                std::move(all_mesh_subsets), vec_n_components,
                NumLib::ComponentOrder::BY_LOCATION);
        assert(_local_to_global_index_map);
    }
    else
    {
        // For displacement equation.
        std::vector<MeshLib::MeshSubsets> all_mesh_subsets;
        std::generate_n(
            std::back_inserter(all_mesh_subsets),
            getProcessVariables(_mechanics_related_process_id)[0]
                .get()
                .getNumberOfComponents(),
            [&]() {
                return MeshLib::MeshSubsets{_mesh_subset_all_nodes.get()};
            });

        std::vector<int> const vec_n_components{DisplacementDim};
        _local_to_global_index_map =
            std::make_unique<NumLib::LocalToGlobalIndexMap>(
                std::move(all_mesh_subsets), vec_n_components,
                NumLib::ComponentOrder::BY_LOCATION);

        // For phase field equation or the heat conduction.
        _sparsity_pattern_with_single_component =
            NumLib::computeSparsityPattern(
                *_local_to_global_index_map_single_component, _mesh);
    }
}

template <int DisplacementDim>
void ThermoMechanicalPhaseFieldProcess<DisplacementDim>::initializeConcreteProcess(
    NumLib::LocalToGlobalIndexMap const& dof_table,
    MeshLib::Mesh const& mesh,
    unsigned const integration_order)
{
    ProcessLib::SmallDeformation::createLocalAssemblers<
        DisplacementDim, ThermoMechanicalPhaseFieldLocalAssembler>(
        mesh.getElements(), dof_table, _local_assemblers,
        mesh.isAxiallySymmetric(), integration_order, _process_data,
        _mechanics_related_process_id, _phase_field_process_id,
        _heat_conduction_process_id);

    Base::_secondary_variables.addSecondaryVariable(
        "sigma_xx",
        makeExtrapolator(
            1, getExtrapolator(), _local_assemblers,
            &ThermoMechanicalPhaseFieldLocalAssemblerInterface::getIntPtSigmaXX));

    Base::_secondary_variables.addSecondaryVariable(
        "sigma_yy",
        makeExtrapolator(
            1, getExtrapolator(), _local_assemblers,
            &ThermoMechanicalPhaseFieldLocalAssemblerInterface::getIntPtSigmaYY));

    Base::_secondary_variables.addSecondaryVariable(
        "sigma_zz",
        makeExtrapolator(
            1, getExtrapolator(), _local_assemblers,
            &ThermoMechanicalPhaseFieldLocalAssemblerInterface::getIntPtSigmaZZ));

    Base::_secondary_variables.addSecondaryVariable(
        "sigma_xy",
        makeExtrapolator(
            1, getExtrapolator(), _local_assemblers,
            &ThermoMechanicalPhaseFieldLocalAssemblerInterface::getIntPtSigmaXY));

    if (DisplacementDim == 3)
    {
        Base::_secondary_variables.addSecondaryVariable(
            "sigma_xz",
            makeExtrapolator(
                1, getExtrapolator(), _local_assemblers,
                &ThermoMechanicalPhaseFieldLocalAssemblerInterface::getIntPtSigmaXZ));

        Base::_secondary_variables.addSecondaryVariable(
            "sigma_yz",
            makeExtrapolator(
                1, getExtrapolator(), _local_assemblers,
                &ThermoMechanicalPhaseFieldLocalAssemblerInterface::getIntPtSigmaYZ));
    }

    Base::_secondary_variables.addSecondaryVariable(
        "epsilon_xx",
        makeExtrapolator(
            1, getExtrapolator(), _local_assemblers,
            &ThermoMechanicalPhaseFieldLocalAssemblerInterface::getIntPtEpsilonXX));

    Base::_secondary_variables.addSecondaryVariable(
        "epsilon_yy",
        makeExtrapolator(
            1, getExtrapolator(), _local_assemblers,
            &ThermoMechanicalPhaseFieldLocalAssemblerInterface::getIntPtEpsilonYY));

    Base::_secondary_variables.addSecondaryVariable(
        "epsilon_zz",
        makeExtrapolator(
            1, getExtrapolator(), _local_assemblers,
            &ThermoMechanicalPhaseFieldLocalAssemblerInterface::getIntPtEpsilonZZ));

    Base::_secondary_variables.addSecondaryVariable(
        "epsilon_xy",
        makeExtrapolator(
            1, getExtrapolator(), _local_assemblers,
            &ThermoMechanicalPhaseFieldLocalAssemblerInterface::getIntPtEpsilonXY));
    if (DisplacementDim == 3)
    {
        Base::_secondary_variables.addSecondaryVariable(
            "epsilon_yz",
            makeExtrapolator(
                1, getExtrapolator(), _local_assemblers,
                &ThermoMechanicalPhaseFieldLocalAssemblerInterface::getIntPtEpsilonYZ));

        Base::_secondary_variables.addSecondaryVariable(
            "epsilon_xz",
            makeExtrapolator(
                1, getExtrapolator(), _local_assemblers,
                &ThermoMechanicalPhaseFieldLocalAssemblerInterface::getIntPtEpsilonXZ));
    }

    Base::_secondary_variables.addSecondaryVariable(
        "heat_flux",
        makeExtrapolator(
            mesh.getDimension(), getExtrapolator(), _local_assemblers,
            &ThermoMechanicalPhaseFieldLocalAssemblerInterface::getIntPtHeatFlux));
}

template <int DisplacementDim>
void ThermoMechanicalPhaseFieldProcess<DisplacementDim>::initializeBoundaryConditions()
{
    if (_use_monolithic_scheme)
    {
        const int monolithic_coupled_processes_id = 0;
        initializeProcessBoundaryConditionsAndSourceTerms(
            *_local_to_global_index_map, monolithic_coupled_processes_id);
        return;
    }

    // Staggered scheme:
    // for the equations of temperature-deformation.
    initializeProcessBoundaryConditionsAndSourceTerms(
        getDOFTableByProcessID(_mechanics_related_process_id),
        _mechanics_related_process_id);
    // for the phase field
    initializeProcessBoundaryConditionsAndSourceTerms(
        getDOFTableByProcessID(_phase_field_process_id),
        _phase_field_process_id);
    // for heat conduction
    initializeProcessBoundaryConditionsAndSourceTerms(
        getDOFTableByProcessID(_heat_conduction_process_id),
        _heat_conduction_process_id);
}

template <int DisplacementDim>
void ThermoMechanicalPhaseFieldProcess<DisplacementDim>::assembleConcreteProcess(
    const double t, GlobalVector const& x, GlobalMatrix& M, GlobalMatrix& K,
    GlobalVector& b)
{
    DBUG("Assemble the equations for ThermoMechanicalPhaseFieldProcess.");

    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_table = {std::ref(*_local_to_global_index_map)};
    // Call global assembler for each local assembly item.
    GlobalExecutor::executeMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assemble, _local_assemblers,
        dof_table, t, x, M, K, b, _coupled_solutions);
}

template <int DisplacementDim>
void ThermoMechanicalPhaseFieldProcess<DisplacementDim>::assembleWithJacobianConcreteProcess(
    const double t, GlobalVector const& x, GlobalVector const& xdot,
    const double dxdot_dx, const double dx_dx, GlobalMatrix& M, GlobalMatrix& K,
    GlobalVector& b, GlobalMatrix& Jac)
{
    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_tables;
    // For the monolithic scheme
    if (_use_monolithic_scheme)
    {
        DBUG("AssembleJacobian ThermoMechanicalPhaseFieldProcess for the monolithic scheme.");
        dof_tables.emplace_back(*_local_to_global_index_map);
    }
    else
    {
        // For the staggered scheme
        if (_coupled_solutions->process_id == _mechanics_related_process_id)
        {
            DBUG(
                "Assemble the Jacobian equations of "
                "temperature-deformation in "
                "ThermoMechanicalPhaseFieldProcess for the staggered scheme.");
        }

        if (_coupled_solutions->process_id == _phase_field_process_id)
        {
            DBUG(
                "Assemble the Jacobian equations of phase field in "
                "ThermoMechanicalPhaseFieldProcess for the staggered scheme.");
        }
        else
        {
            DBUG(
                "Assemble the Jacobian equations of "
                "heat conduction in "
                "ThermoMechanicalPhaseFieldProcess for the staggered scheme.");
        }
        dof_tables.emplace_back(
            getDOFTableByProcessID(_heat_conduction_process_id));
        dof_tables.emplace_back(
            getDOFTableByProcessID(_mechanics_related_process_id));
        dof_tables.emplace_back(
            getDOFTableByProcessID(_phase_field_process_id));
    }

    GlobalExecutor::executeMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assembleWithJacobian,
        _local_assemblers, dof_tables, t, x, xdot, dxdot_dx, dx_dx, M, K, b,
        Jac, _coupled_solutions);
}

template <int DisplacementDim>
void ThermoMechanicalPhaseFieldProcess<DisplacementDim>::preTimestepConcreteProcess(
    GlobalVector const& x, double const t, double const dt,
    const int process_id)
{
    DBUG("PreTimestep ThermoMechanicalPhaseFieldProcess.");

    _process_data.dt = dt;
    _process_data.t = t;

    if (process_id != _mechanics_related_process_id)
    {
        return;
    }
    GlobalExecutor::executeMemberOnDereferenced(
        &ThermoMechanicalPhaseFieldLocalAssemblerInterface::preTimestep, _local_assemblers,
        getDOFTable(process_id), x, t, dt);
}

template <int DisplacementDim>
void ThermoMechanicalPhaseFieldProcess<DisplacementDim>::postTimestepConcreteProcess(
    GlobalVector const& x, int const process_id)
{
    DBUG("PostTimestep ThermoMechanicalPhaseFieldProcess.");

    GlobalExecutor::executeMemberOnDereferenced(
        &ThermoMechanicalPhaseFieldLocalAssemblerInterface::postTimestep, _local_assemblers,
        getDOFTable(process_id), x);
}

template <int DisplacementDim>
void ThermoMechanicalPhaseFieldProcess<DisplacementDim>::postNonLinearSolverConcreteProcess(
    GlobalVector const& x, const double t, const int process_id)
{
    if (process_id != _mechanics_related_process_id)
    {
        return;
    }

    DBUG("PostNonLinearSolver ThermoMechanicalPhaseFieldProcess.");
    // Calculate strain, stress or other internal variables of mechanics.
    GlobalExecutor::executeMemberOnDereferenced(
        &LocalAssemblerInterface::postNonLinearSolver, _local_assemblers,
        getDOFTable(process_id), x, t, _use_monolithic_scheme);
}

}  // namespace ThermoMechanicalPhaseField
}  // namespace ProcessLib
