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
    // Disable nodal forces for monolithic scheme.
    if (!_use_monolithic_scheme)
    {
        _nodal_forces = MeshLib::getOrCreateMeshProperty<double>(
            mesh, "NodalForces", MeshLib::MeshItemType::Node, DisplacementDim);
    }
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
    if (_use_monolithic_scheme || process_id == 0)
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
    if (_use_monolithic_scheme || process_id == 0)
    {
        return *_local_to_global_index_map;
    }

    // For the equation of phasefield
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
        const int process_id = 0;
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
        "sigma",
        makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                         &LocalAssemblerInterface::getIntPtSigma));

    Base::_secondary_variables.addSecondaryVariable(
        "epsilon",
        makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                         &LocalAssemblerInterface::getIntPtEpsilon));
}

template <int DisplacementDim>
void PhaseFieldProcess<DisplacementDim>::initializeBoundaryConditions()
{
    if (_use_monolithic_scheme)
    {
        const int process_id_of_phasefieldmechanics = 0;
        initializeProcessBoundaryConditionsAndSourceTerms(
            *_local_to_global_index_map, process_id_of_phasefieldmechanics);
        return;
    }

    // Staggered scheme:
    // for the equations of deformation.
    const int mechanical_process_id = 0;
    initializeProcessBoundaryConditionsAndSourceTerms(
        *_local_to_global_index_map, mechanical_process_id);
    // for the phase field
    const int phasefield_process_id = 1;
    initializeProcessBoundaryConditionsAndSourceTerms(
        *_local_to_global_index_map_single_component, phasefield_process_id);
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
        if (_coupled_solutions->process_id == 1)
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
        if (_coupled_solutions->process_id == 1)
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
        dof_tables.emplace_back(*_local_to_global_index_map);
        dof_tables.emplace_back(*_local_to_global_index_map_single_component);
    }
    // Call global assembler for each local assembly item.

    GlobalExecutor::executeMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assembleWithJacobian,
        _local_assemblers, dof_tables, t, x, xdot, dxdot_dx, dx_dx, M, K, b,
        Jac, _coupled_solutions);

    if (!_use_monolithic_scheme && (_coupled_solutions->process_id == 0))
    {
        b.copyValues(*_nodal_forces);
        std::transform(_nodal_forces->begin(), _nodal_forces->end(),
                       _nodal_forces->begin(), [](double val) { return -val; });
    }
}

template <int DisplacementDim>
void PhaseFieldProcess<DisplacementDim>::preTimestepConcreteProcess(
    GlobalVector const& x, double const t, double const dt,
    const int process_id)
{
    DBUG("PreTimestep PhaseFieldProcess %d.", process_id);

    _process_data.dt = dt;
    _process_data.t = t;
    _process_data.injected_volume = _process_data.t;

    GlobalExecutor::executeMemberOnDereferenced(
        &LocalAssemblerInterface::preTimestep, _local_assemblers,
        getDOFTable(process_id), x, t, dt);
}

template <int DisplacementDim>
void PhaseFieldProcess<DisplacementDim>::postTimestepConcreteProcess(
    GlobalVector const& x, int const process_id)
{
    if (isPhaseFieldProcess(process_id))
    {
        DBUG("PostTimestep PhaseFieldProcess.");

        _process_data.elastic_energy = 0.0;
        _process_data.surface_energy = 0.0;
        _process_data.pressure_work = 0.0;

        std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
            dof_tables;

        dof_tables.emplace_back(*_local_to_global_index_map);
        dof_tables.emplace_back(*_local_to_global_index_map_single_component);

        GlobalExecutor::executeMemberOnDereferenced(
            &LocalAssemblerInterface::computeEnergy, _local_assemblers,
            dof_tables, x, _process_data.t, _process_data.elastic_energy,
            _process_data.surface_energy, _process_data.pressure_work,
            _use_monolithic_scheme, _coupled_solutions);

        INFO("Elastic energy: %g Surface energy: %g Pressure work: %g ",
             _process_data.elastic_energy, _process_data.surface_energy,
             _process_data.pressure_work);
    }
}

template <int DisplacementDim>
void PhaseFieldProcess<DisplacementDim>::postNonLinearSolverConcreteProcess(
    GlobalVector const& x, const double t, const int process_id)
{
    if (_use_monolithic_scheme)
        return;

    _process_data.crack_volume = 0.0;

    if (!isPhaseFieldProcess(process_id))
    {
        std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
            dof_tables;

        dof_tables.emplace_back(*_local_to_global_index_map);
        dof_tables.emplace_back(*_local_to_global_index_map_single_component);

        DBUG("PostNonLinearSolver crack volume computation.");

        GlobalExecutor::executeMemberOnDereferenced(
            &LocalAssemblerInterface::computeCrackIntegral, _local_assemblers,
            dof_tables, x, t, _process_data.crack_volume,
            _use_monolithic_scheme, _coupled_solutions);

        INFO("Integral of crack: %g", _process_data.crack_volume);

        if (_process_data.propagating_crack)
        {
            _process_data.pressure_old = _process_data.pressure;
            _process_data.pressure =
                _process_data.injected_volume / _process_data.crack_volume;
            _process_data.pressure_error =
                abs(_process_data.pressure_old - _process_data.pressure) /
                _process_data.pressure;
            INFO("Internal pressure: %g and Pressure error: %.4e",
                 _process_data.pressure, _process_data.pressure_error);
            auto& u = _coupled_solutions->coupled_xs[0].get();
            MathLib::LinAlg::scale(const_cast<GlobalVector&>(u),
                                   _process_data.pressure);
        }
    }
    else
    {
        if (_process_data.propagating_crack)
        {
            auto& u = _coupled_solutions->coupled_xs[0].get();
            MathLib::LinAlg::scale(const_cast<GlobalVector&>(u),
                                   1 / _process_data.pressure);
        }
    }
}
}  // namespace PhaseField
}  // namespace ProcessLib
