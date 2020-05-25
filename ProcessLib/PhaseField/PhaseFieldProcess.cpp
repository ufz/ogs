/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "PhaseFieldProcess.h"

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
    std::string name,
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    unsigned const integration_order,
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
        process_variables,
    PhaseFieldProcessData<DisplacementDim>&& process_data,
    SecondaryVariableCollection&& secondary_variables,
    bool const use_monolithic_scheme)
    : Process(std::move(name), mesh, std::move(jacobian_assembler), parameters,
              integration_order, std::move(process_variables),
              std::move(secondary_variables), use_monolithic_scheme),
      process_data_(std::move(process_data))
{
    if (use_monolithic_scheme)
    {
        OGS_FATAL(
            "Monolithic scheme is not implemented for the PhaseField process.");
    }

    nodal_forces_ = MeshLib::getOrCreateMeshProperty<double>(
        mesh, "NodalForces", MeshLib::MeshItemType::Node, DisplacementDim);
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
    // For the M process (deformation) in the staggered scheme.
    if (process_id == 0)
    {
        auto const& l = *local_to_global_index_map_;
        return {l.dofSizeWithoutGhosts(), l.dofSizeWithoutGhosts(),
                &l.getGhostIndices(), &this->sparsity_pattern_};
    }

    // For staggered scheme and phase field process.
    auto const& l = *local_to_global_index_map_single_component_;
    return {l.dofSizeWithoutGhosts(), l.dofSizeWithoutGhosts(),
            &l.getGhostIndices(), &sparsity_pattern_with_single_component_};
}

template <int DisplacementDim>
NumLib::LocalToGlobalIndexMap const&
PhaseFieldProcess<DisplacementDim>::getDOFTable(const int process_id) const
{
    // For the M process (deformation) in the staggered scheme.
    if (process_id == 0)
    {
        return *local_to_global_index_map_;
    }

    // For the equation of phasefield
    return *local_to_global_index_map_single_component_;
}

template <int DisplacementDim>
void PhaseFieldProcess<DisplacementDim>::constructDofTable()
{
    // For displacement equation.
    const int mechanics_process_id = 0;
    constructDofTableOfSpecifiedProsessStaggerdScheme(mechanics_process_id);

    // TODO move the two data members somewhere else.
    // for extrapolation of secondary variables of stress or strain
    std::vector<MeshLib::MeshSubset> all_mesh_subsets_single_component{
        *mesh_subset_all_nodes_};
    local_to_global_index_map_single_component_ =
        std::make_unique<NumLib::LocalToGlobalIndexMap>(
            std::move(all_mesh_subsets_single_component),
            // by location order is needed for output
            NumLib::ComponentOrder::BY_LOCATION);

    assert(local_to_global_index_map_single_component_);

    // For phase field equation.
    sparsity_pattern_with_single_component_ = NumLib::computeSparsityPattern(
        *local_to_global_index_map_single_component_, mesh_);
}

template <int DisplacementDim>
void PhaseFieldProcess<DisplacementDim>::initializeConcreteProcess(
    NumLib::LocalToGlobalIndexMap const& dof_table,
    MeshLib::Mesh const& mesh,
    unsigned const integration_order)
{
    ProcessLib::SmallDeformation::createLocalAssemblers<
        DisplacementDim, PhaseFieldLocalAssembler>(
        mesh.getElements(), dof_table, local_assemblers_,
        mesh.isAxiallySymmetric(), integration_order, process_data_);

    secondary_variables_.addSecondaryVariable(
        "sigma",
        makeExtrapolator(MathLib::KelvinVector::KelvinVectorType<
                             DisplacementDim>::RowsAtCompileTime,
                         getExtrapolator(), local_assemblers_,
                         &LocalAssemblerInterface::getIntPtSigma));

    secondary_variables_.addSecondaryVariable(
        "epsilon",
        makeExtrapolator(MathLib::KelvinVector::KelvinVectorType<
                             DisplacementDim>::RowsAtCompileTime,
                         getExtrapolator(), local_assemblers_,
                         &LocalAssemblerInterface::getIntPtEpsilon));

    // Initialize local assemblers after all variables have been set.
    GlobalExecutor::executeMemberOnDereferenced(
        &LocalAssemblerInterface::initialize, local_assemblers_,
        *local_to_global_index_map_);
}

template <int DisplacementDim>
void PhaseFieldProcess<DisplacementDim>::initializeBoundaryConditions()
{
    // Staggered scheme:
    // for the equations of deformation.
    const int mechanical_process_id = 0;
    initializeProcessBoundaryConditionsAndSourceTerms(
        *local_to_global_index_map_, mechanical_process_id);
    // for the phase field
    const int phasefield_process_id = 1;
    initializeProcessBoundaryConditionsAndSourceTerms(
        *local_to_global_index_map_single_component_, phasefield_process_id);
}

template <int DisplacementDim>
void PhaseFieldProcess<DisplacementDim>::assembleConcreteProcess(
    const double t, double const dt, std::vector<GlobalVector*> const& x,
    std::vector<GlobalVector*> const& xdot, int const process_id,
    GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b)
{
    DBUG("Assemble PhaseFieldProcess.");

    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_tables;

    // For the staggered scheme
    if (process_id == 1)
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
    dof_tables.emplace_back(*local_to_global_index_map_single_component_);
    dof_tables.emplace_back(*local_to_global_index_map_);

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    // Call global assembler for each local assembly item.
    GlobalExecutor::executeSelectedMemberDereferenced(
        global_assembler_, &VectorMatrixAssembler::assemble, local_assemblers_,
        pv.getActiveElementIDs(), dof_tables, t, dt, x, xdot, process_id, M, K,
        b, coupled_solutions_);
}

template <int DisplacementDim>
void PhaseFieldProcess<DisplacementDim>::assembleWithJacobianConcreteProcess(
    const double t, double const dt, std::vector<GlobalVector*> const& x,
    GlobalVector const& xdot, const double dxdot_dx, const double dx_dx,
    int const process_id, GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b,
    GlobalMatrix& Jac)
{
    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_tables;

    // For the staggered scheme
    if (process_id == 1)
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
    dof_tables.emplace_back(*local_to_global_index_map_);
    dof_tables.emplace_back(*local_to_global_index_map_single_component_);

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    // Call global assembler for each local assembly item.
    GlobalExecutor::executeSelectedMemberDereferenced(
        global_assembler_, &VectorMatrixAssembler::assembleWithJacobian,
        local_assemblers_, pv.getActiveElementIDs(), dof_tables, t, dt, x, xdot,
        dxdot_dx, dx_dx, process_id, M, K, b, Jac, coupled_solutions_);

    if (process_id == 0)
    {
        b.copyValues(*nodal_forces_);
        std::transform(nodal_forces_->begin(), nodal_forces_->end(),
                       nodal_forces_->begin(), [](double val) { return -val; });
    }
}

template <int DisplacementDim>
void PhaseFieldProcess<DisplacementDim>::preTimestepConcreteProcess(
    std::vector<GlobalVector*> const& x, double const t, double const dt,
    const int process_id)
{
    DBUG("PreTimestep PhaseFieldProcess {:d}.", process_id);

    process_data_.injected_volume = t;

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    GlobalExecutor::executeSelectedMemberOnDereferenced(
        &LocalAssemblerInterface::preTimestep, local_assemblers_,
        pv.getActiveElementIDs(), getDOFTable(process_id), *x[process_id], t,
        dt);
}

template <int DisplacementDim>
void PhaseFieldProcess<DisplacementDim>::postTimestepConcreteProcess(
    std::vector<GlobalVector*> const& x, const double t,
    const double /*delta_t*/, int const process_id)
{
    if (isPhaseFieldProcess(process_id))
    {
        DBUG("PostTimestep PhaseFieldProcess.");

        process_data_.elastic_energy = 0.0;
        process_data_.surface_energy = 0.0;
        process_data_.pressure_work = 0.0;

        std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
            dof_tables;

        dof_tables.emplace_back(*local_to_global_index_map_);
        dof_tables.emplace_back(*local_to_global_index_map_single_component_);

        ProcessLib::ProcessVariable const& pv =
            getProcessVariables(process_id)[0];

        GlobalExecutor::executeSelectedMemberOnDereferenced(
            &LocalAssemblerInterface::computeEnergy, local_assemblers_,
            pv.getActiveElementIDs(), dof_tables, *x[process_id], t,
            process_data_.elastic_energy, process_data_.surface_energy,
            process_data_.pressure_work, coupled_solutions_);

        INFO("Elastic energy: {:g} Surface energy: {:g} Pressure work: {:g} ",
             process_data_.elastic_energy, process_data_.surface_energy,
             process_data_.pressure_work);
    }
}

template <int DisplacementDim>
void PhaseFieldProcess<DisplacementDim>::postNonLinearSolverConcreteProcess(
    GlobalVector const& x, const double t, double const /*dt*/,
    const int process_id)
{
    process_data_.crack_volume = 0.0;

    if (!isPhaseFieldProcess(process_id))
    {
        std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
            dof_tables;

        dof_tables.emplace_back(*local_to_global_index_map_);
        dof_tables.emplace_back(*local_to_global_index_map_single_component_);

        DBUG("PostNonLinearSolver crack volume computation.");

        ProcessLib::ProcessVariable const& pv =
            getProcessVariables(process_id)[0];
        GlobalExecutor::executeSelectedMemberOnDereferenced(
            &LocalAssemblerInterface::computeCrackIntegral, local_assemblers_,
            pv.getActiveElementIDs(), dof_tables, x, t,
            process_data_.crack_volume, coupled_solutions_);

        INFO("Integral of crack: {:g}", process_data_.crack_volume);

        if (process_data_.propagating_crack)
        {
            process_data_.pressure_old = process_data_.pressure;
            process_data_.pressure =
                process_data_.injected_volume / process_data_.crack_volume;
            process_data_.pressure_error =
                std::fabs(process_data_.pressure_old - process_data_.pressure) /
                process_data_.pressure;
            INFO("Internal pressure: {:g} and Pressure error: {:.4e}",
                 process_data_.pressure, process_data_.pressure_error);
            auto& u = *coupled_solutions_->coupled_xs[0];
            MathLib::LinAlg::scale(const_cast<GlobalVector&>(u),
                                   process_data_.pressure);
        }
    }
    else
    {
        if (process_data_.propagating_crack)
        {
            auto& u = *coupled_solutions_->coupled_xs[0];
            MathLib::LinAlg::scale(const_cast<GlobalVector&>(u),
                                   1 / process_data_.pressure);
        }
    }
}

template <int DisplacementDim>
bool PhaseFieldProcess<DisplacementDim>::isPhaseFieldProcess(
    int const process_id) const
{
    return process_id == 1;
}

template class PhaseFieldProcess<2>;
template class PhaseFieldProcess<3>;

}  // namespace PhaseField
}  // namespace ProcessLib
