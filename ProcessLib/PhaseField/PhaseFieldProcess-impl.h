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

#include "NumLib/DOF/ComputeSparsityPattern.h"

#include "ProcessLib/Process.h"
#include "ProcessLib/SmallDeformation/CreateLocalAssemblers.h"

#include "PhaseFieldFEM.h"

namespace ProcessLib
{
namespace PhaseField
{
template <int DisplacementDim>
PhaseFieldProcess<DisplacementDim>::PhaseFieldProcess(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order,
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
        process_variables,
    PhaseFieldProcessData<DisplacementDim>&& process_data,
    SecondaryVariableCollection&& secondary_variables,
    NumLib::NamedFunctionCaller&& named_function_caller,
    bool const use_monolithic_scheme)
    : Process(mesh, std::move(jacobian_assembler), parameters,
              integration_order, std::move(process_variables),
              std::move(secondary_variables), std::move(named_function_caller),
              use_monolithic_scheme),
      _process_data(std::move(process_data))
{
}

template <int DisplacementDim>
bool PhaseFieldProcess<DisplacementDim>::isLinear() const
{
    return false;
}

template <int DisplacementDim>
MathLib::MatrixSpecifications
PhaseFieldProcess<DisplacementDim>::getMatrixSpecifications(
    const int process_id) const
{
    // For the monolithic scheme or the M process (deformation) in the staggered
    // scheme.
    if (_use_monolithic_scheme || process_id == 1)
    {
        auto const& l = *_local_to_global_index_map;
        return {l.dofSizeWithoutGhosts(), l.dofSizeWithoutGhosts(),
                &l.getGhostIndices(), &this->_sparsity_pattern};
    }

    // For staggered scheme and phase field process.
    auto const& l = *_local_to_global_index_map_single_component;
    return {l.dofSizeWithoutGhosts(), l.dofSizeWithoutGhosts(),
            &l.getGhostIndices(), &_sparsity_pattern_with_single_component};
}

template <int DisplacementDim>
NumLib::LocalToGlobalIndexMap const&
PhaseFieldProcess<DisplacementDim>::getDOFTable(const int process_id) const
{
    // If monolithic scheme is used or the equation of deformation is solved in
    // the staggered scheme.
    if (_use_monolithic_scheme || process_id == 1)
    {
        return *_local_to_global_index_map;
    }

    // For the equation of pressure
    return *_local_to_global_index_map_single_component;
}

template <int DisplacementDim>
void PhaseFieldProcess<DisplacementDim>::constructDofTable()
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
        const int process_id = 0;  // Only one process in the monolithic scheme.
        std::vector<MeshLib::MeshSubsets> all_mesh_subsets;
        std::generate_n(
            std::back_inserter(all_mesh_subsets),
            getProcessVariables(process_id)[1].get().getNumberOfComponents() +
                1,
            [&]() {
                return MeshLib::MeshSubsets{_mesh_subset_all_nodes.get()};
            });

        std::vector<int> const vec_n_components{1, DisplacementDim};
        _local_to_global_index_map =
            std::make_unique<NumLib::LocalToGlobalIndexMap>(
                std::move(all_mesh_subsets), vec_n_components,
                NumLib::ComponentOrder::BY_LOCATION);
        assert(_local_to_global_index_map);
    }
    else
    {
        // For displacement equation.
        const int process_id = 1;
        std::vector<MeshLib::MeshSubsets> all_mesh_subsets;
        std::generate_n(
            std::back_inserter(all_mesh_subsets),
            getProcessVariables(process_id)[0].get().getNumberOfComponents(),
            [&]() {
                return MeshLib::MeshSubsets{_mesh_subset_all_nodes.get()};
            });

        std::vector<int> const vec_n_components{DisplacementDim};
        _local_to_global_index_map =
            std::make_unique<NumLib::LocalToGlobalIndexMap>(
                std::move(all_mesh_subsets), vec_n_components,
                NumLib::ComponentOrder::BY_LOCATION);

        // For phase field equation.
        _sparsity_pattern_with_single_component =
            NumLib::computeSparsityPattern(
                *_local_to_global_index_map_single_component, _mesh);
    }
}

template <int DisplacementDim>
void PhaseFieldProcess<DisplacementDim>::initializeConcreteProcess(
    NumLib::LocalToGlobalIndexMap const& dof_table,
    MeshLib::Mesh const& mesh,
    unsigned const integration_order)
{
    ProcessLib::SmallDeformation::createLocalAssemblers<
        DisplacementDim, PhaseFieldLocalAssembler>(
        mesh.getElements(), dof_table, _local_assemblers,
        mesh.isAxiallySymmetric(), integration_order, _process_data);

    Base::_secondary_variables.addSecondaryVariable(
        "sigma_xx",
        makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                         &LocalAssemblerInterface::getIntPtSigmaXX));

    Base::_secondary_variables.addSecondaryVariable(
        "sigma_yy",
        makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                         &LocalAssemblerInterface::getIntPtSigmaYY));

    Base::_secondary_variables.addSecondaryVariable(
        "sigma_zz",
        makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                         &LocalAssemblerInterface::getIntPtSigmaZZ));

    Base::_secondary_variables.addSecondaryVariable(
        "sigma_xy",
        makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                         &LocalAssemblerInterface::getIntPtSigmaXY));

    if (DisplacementDim == 3)
    {
        Base::_secondary_variables.addSecondaryVariable(
            "sigma_xz",
            makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                             &LocalAssemblerInterface::getIntPtSigmaXZ));

        Base::_secondary_variables.addSecondaryVariable(
            "sigma_yz",
            makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                             &LocalAssemblerInterface::getIntPtSigmaYZ));
    }

    Base::_secondary_variables.addSecondaryVariable(
        "epsilon_xx",
        makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                         &LocalAssemblerInterface::getIntPtEpsilonXX));

    Base::_secondary_variables.addSecondaryVariable(
        "epsilon_yy",
        makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                         &LocalAssemblerInterface::getIntPtEpsilonYY));

    Base::_secondary_variables.addSecondaryVariable(
        "epsilon_zz",
        makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                         &LocalAssemblerInterface::getIntPtEpsilonZZ));

    Base::_secondary_variables.addSecondaryVariable(
        "epsilon_xy",
        makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                         &LocalAssemblerInterface::getIntPtEpsilonXY));

    if (DisplacementDim == 3)
    {
        Base::_secondary_variables.addSecondaryVariable(
            "epsilon_yz",
            makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                             &LocalAssemblerInterface::getIntPtEpsilonYZ));

        Base::_secondary_variables.addSecondaryVariable(
            "epsilon_xz",
            makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                             &LocalAssemblerInterface::getIntPtEpsilonXZ));
    }
}

template <int DisplacementDim>
void PhaseFieldProcess<DisplacementDim>::initializeBoundaryConditions()
{
    if (_use_monolithic_scheme)
    {
        const int process_id_of_pf = 0;
        initializeProcessBoundaryConditionsAndSourceTerms(
            *_local_to_global_index_map, process_id_of_pf);
        return;
    }

    // Staggered scheme:
    // for the phase field
    const int process_id_of_pf = 0;
    initializeProcessBoundaryConditionsAndSourceTerms(
        *_local_to_global_index_map_single_component, process_id_of_pf);

    // for the equations of deformation.
    const int process_id_of_u = 1;
    initializeProcessBoundaryConditionsAndSourceTerms(
        *_local_to_global_index_map, process_id_of_u);
}

template <int DisplacementDim>
void PhaseFieldProcess<DisplacementDim>::assembleConcreteProcess(
    const double t, GlobalVector const& x, GlobalMatrix& M, GlobalMatrix& K,
    GlobalVector& b)
{
    DBUG("Assemble PhaseFieldProcess.");

    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_tables;

    // For the monolithic scheme
    if (_use_monolithic_scheme)
    {
        DBUG("Assemble PhaseFieldProcess for the monolithic scheme.");
        dof_tables.emplace_back(*_local_to_global_index_map);
    }
    else
    {
        // For the staggered scheme
        if (_coupled_solutions->process_id == 0)
        {
            DBUG(
                "Assemble the equations of phase field in "
                "PhaseFieldProcess for the staggered scheme.");
        }
        else
        {
            DBUG(
                "Assemble the equations of deformation in "
                "PhaseFieldProcess for the staggered scheme.");
        }
        dof_tables.emplace_back(*_local_to_global_index_map_single_component);
        dof_tables.emplace_back(*_local_to_global_index_map);
    }

    // Call global assembler for each local assembly item.
    GlobalExecutor::executeMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assemble, _local_assemblers,
        dof_tables, t, x, M, K, b, _coupled_solutions);
}

template <int DisplacementDim>
void PhaseFieldProcess<DisplacementDim>::assembleWithJacobianConcreteProcess(
    const double t, GlobalVector const& x, GlobalVector const& xdot,
    const double dxdot_dx, const double dx_dx, GlobalMatrix& M, GlobalMatrix& K,
    GlobalVector& b, GlobalMatrix& Jac)
{
    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_tables;

    // For the monolithic scheme
    if (_use_monolithic_scheme)
    {
        DBUG("AssembleJacobian PhaseFieldProcess for the monolithic scheme.");
        dof_tables.emplace_back(*_local_to_global_index_map);
    }
    else
    {
        // For the staggered scheme
        if (_coupled_solutions->process_id == 0)
        {
            DBUG(
                "Assemble the Jacobian equations of phase field in "
                "PhaseFieldProcess for the staggered scheme.");
        }
        else
        {
            DBUG(
                "Assemble the Jacobian equations of deformation in "
                "PhaseFieldProcess for the staggered scheme.");
        }
        dof_tables.emplace_back(*_local_to_global_index_map_single_component);
        dof_tables.emplace_back(*_local_to_global_index_map);
    }
    // Call global assembler for each local assembly item.
    GlobalExecutor::executeMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assembleWithJacobian,
        _local_assemblers, dof_tables, t, x, xdot, dxdot_dx, dx_dx, M, K, b,
        Jac, _coupled_solutions);
}

template <int DisplacementDim>
void PhaseFieldProcess<DisplacementDim>::preTimestepConcreteProcess(
    GlobalVector const& x, double const t, double const dt,
    const int process_id)
{
    DBUG("PreTimestep PhaseFieldProcess.");

    _process_data.dt = dt;
    _process_data.t = t;

    GlobalExecutor::executeMemberOnDereferenced(
        &LocalAssemblerInterface::preTimestep, _local_assemblers,
        getDOFTable(process_id), x, t, dt);
}

template <int DisplacementDim>
void PhaseFieldProcess<DisplacementDim>::postTimestepConcreteProcess(
    GlobalVector const& x, int const process_id)
{
    DBUG("PostTimestep PhaseFieldProcess.");

    GlobalExecutor::executeMemberOnDereferenced(
        &LocalAssemblerInterface::postTimestep, _local_assemblers,
        getDOFTable(process_id), x);
}

template <int DisplacementDim>
void PhaseFieldProcess<DisplacementDim>::postNonLinearSolverConcreteProcess(
    GlobalVector const& x, const double t, const int process_id)
{
    if (!_use_monolithic_scheme && process_id == 0)
    {
        return;
    }

    DBUG("PostNonLinearSolver PhaseFieldProcess.");
    // Calculate strain, stress or other internal variables of mechanics.
    GlobalExecutor::executeMemberOnDereferenced(
        &LocalAssemblerInterface::postNonLinearSolver, _local_assemblers,
        getDOFTable(process_id), x, t, _use_monolithic_scheme);
}
}  // namespace PhaseField
}  // namespace ProcessLib
