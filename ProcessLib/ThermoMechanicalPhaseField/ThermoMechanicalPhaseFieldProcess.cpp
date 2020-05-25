/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ThermoMechanicalPhaseFieldFEM.h"
#include "ThermoMechanicalPhaseFieldProcess.h"
#include "ThermoMechanicalPhaseFieldProcessData.h"

#include <cassert>

#include "NumLib/DOF/ComputeSparsityPattern.h"
#include "ProcessLib/Process.h"
#include "ProcessLib/SmallDeformation/CreateLocalAssemblers.h"

namespace ProcessLib
{
namespace ThermoMechanicalPhaseField
{
template <int DisplacementDim>
ThermoMechanicalPhaseFieldProcess<DisplacementDim>::
    ThermoMechanicalPhaseFieldProcess(
        std::string name,
        MeshLib::Mesh& mesh,
        std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&&
            jacobian_assembler,
        std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const&
            parameters,
        unsigned const integration_order,
        std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
            process_variables,
        ThermoMechanicalPhaseFieldProcessData<DisplacementDim>&& process_data,
        SecondaryVariableCollection&& secondary_variables,
        int const mechanics_related_process_id,
        int const phase_field_process_id,
        int const heat_conduction_process_id)
    : Process(std::move(name), mesh, std::move(jacobian_assembler), parameters,
              integration_order, std::move(process_variables),
              std::move(secondary_variables), false),
      process_data_(std::move(process_data)),
      mechanics_related_process_id_(mechanics_related_process_id),
      phase_field_process_id_(phase_field_process_id),
      heat_conduction_process_id_(heat_conduction_process_id)
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
    if (process_id == mechanics_related_process_id_)
    {
        auto const& l = *local_to_global_index_map_;
        return {l.dofSizeWithoutGhosts(), l.dofSizeWithoutGhosts(),
                &l.getGhostIndices(), &this->sparsity_pattern_};
    }

    // For staggered scheme and phase field process or heat conduction.
    auto const& l = *local_to_global_index_map_single_component_;
    return {l.dofSizeWithoutGhosts(), l.dofSizeWithoutGhosts(),
            &l.getGhostIndices(), &sparsity_pattern_with_single_component_};
}

template <int DisplacementDim>
NumLib::LocalToGlobalIndexMap const&
ThermoMechanicalPhaseFieldProcess<DisplacementDim>::getDOFTable(
    const int process_id) const
{
    if (process_id == mechanics_related_process_id_)
    {
        return *local_to_global_index_map_;
    }

    // For the equation of phasefield or heat conduction.
    return *local_to_global_index_map_single_component_;
}

template <int DisplacementDim>
NumLib::LocalToGlobalIndexMap&
ThermoMechanicalPhaseFieldProcess<DisplacementDim>::getDOFTableByProcessID(
    const int process_id) const
{
    if (process_id == mechanics_related_process_id_)
    {
        return *local_to_global_index_map_;
    }

    // For the equation of phasefield or heat conduction.
    return *local_to_global_index_map_single_component_;
}

template <int DisplacementDim>
void ThermoMechanicalPhaseFieldProcess<DisplacementDim>::constructDofTable()
{
    // For displacement equation.
    constructDofTableOfSpecifiedProsessStaggerdScheme(
        mechanics_related_process_id_);

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

    // For phase field equation or the heat conduction.
    sparsity_pattern_with_single_component_ = NumLib::computeSparsityPattern(
        *local_to_global_index_map_single_component_, mesh_);
}

template <int DisplacementDim>
void ThermoMechanicalPhaseFieldProcess<DisplacementDim>::
    initializeConcreteProcess(NumLib::LocalToGlobalIndexMap const& dof_table,
                              MeshLib::Mesh const& mesh,
                              unsigned const integration_order)
{
    ProcessLib::SmallDeformation::createLocalAssemblers<
        DisplacementDim, ThermoMechanicalPhaseFieldLocalAssembler>(
        mesh.getElements(), dof_table, local_assemblers_,
        mesh.isAxiallySymmetric(), integration_order, process_data_,
        mechanics_related_process_id_, phase_field_process_id_,
        heat_conduction_process_id_);

    secondary_variables_.addSecondaryVariable(
        "sigma",
        makeExtrapolator(
            MathLib::KelvinVector::KelvinVectorType<
                DisplacementDim>::RowsAtCompileTime,
            getExtrapolator(), local_assemblers_,
            &ThermoMechanicalPhaseFieldLocalAssemblerInterface::getIntPtSigma));

    secondary_variables_.addSecondaryVariable(
        "epsilon",
        makeExtrapolator(MathLib::KelvinVector::KelvinVectorType<
                             DisplacementDim>::RowsAtCompileTime,
                         getExtrapolator(), local_assemblers_,
                         &ThermoMechanicalPhaseFieldLocalAssemblerInterface::
                             getIntPtEpsilon));

    secondary_variables_.addSecondaryVariable(
        "heat_flux",
        makeExtrapolator(mesh.getDimension(), getExtrapolator(),
                         local_assemblers_,
                         &ThermoMechanicalPhaseFieldLocalAssemblerInterface::
                             getIntPtHeatFlux));

    // Initialize local assemblers after all variables have been set.
    GlobalExecutor::executeMemberOnDereferenced(
        &LocalAssemblerInterface::initialize, local_assemblers_,
        *local_to_global_index_map_);
}

template <int DisplacementDim>
void ThermoMechanicalPhaseFieldProcess<
    DisplacementDim>::initializeBoundaryConditions()
{
    // Staggered scheme:
    // for the equations of temperature-deformation.
    initializeProcessBoundaryConditionsAndSourceTerms(
        getDOFTableByProcessID(mechanics_related_process_id_),
        mechanics_related_process_id_);
    // for the phase field
    initializeProcessBoundaryConditionsAndSourceTerms(
        getDOFTableByProcessID(phase_field_process_id_),
        phase_field_process_id_);
    // for heat conduction
    initializeProcessBoundaryConditionsAndSourceTerms(
        getDOFTableByProcessID(heat_conduction_process_id_),
        heat_conduction_process_id_);
}

template <int DisplacementDim>
void ThermoMechanicalPhaseFieldProcess<DisplacementDim>::
    assembleConcreteProcess(const double t, double const dt,
                            std::vector<GlobalVector*> const& x,
                            std::vector<GlobalVector*> const& xdot,
                            int const process_id, GlobalMatrix& M,
                            GlobalMatrix& K, GlobalVector& b)
{
    DBUG("Assemble the equations for ThermoMechanicalPhaseFieldProcess.");

    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_table = {std::ref(*local_to_global_index_map_)};
    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    // Call global assembler for each local assembly item.
    GlobalExecutor::executeSelectedMemberDereferenced(
        global_assembler_, &VectorMatrixAssembler::assemble, local_assemblers_,
        pv.getActiveElementIDs(), dof_table, t, dt, x, xdot, process_id, M, K,
        b, coupled_solutions_);
}

template <int DisplacementDim>
void ThermoMechanicalPhaseFieldProcess<DisplacementDim>::
    assembleWithJacobianConcreteProcess(
        const double t, double const dt, std::vector<GlobalVector*> const& x,
        GlobalVector const& xdot, const double dxdot_dx, const double dx_dx,
        int const process_id, GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b,
        GlobalMatrix& Jac)
{
    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_tables;
    // For the staggered scheme
    if (process_id == mechanics_related_process_id_)
    {
        DBUG(
            "Assemble the Jacobian equations of "
            "temperature-deformation in "
            "ThermoMechanicalPhaseFieldProcess for "
            "the staggered scheme.");
    }

    if (process_id == phase_field_process_id_)
    {
        DBUG(
            "Assemble the Jacobian equations of"
            "phase field in "
            "ThermoMechanicalPhaseFieldProcess for "
            "the staggered scheme.");
    }
    else
    {
        DBUG(
            "Assemble the Jacobian equations of "
            "heat conduction in "
            "ThermoMechanicalPhaseFieldProcess for "
            "the staggered scheme.");
    }
    dof_tables.emplace_back(
        getDOFTableByProcessID(heat_conduction_process_id_));
    dof_tables.emplace_back(
        getDOFTableByProcessID(mechanics_related_process_id_));
    dof_tables.emplace_back(getDOFTableByProcessID(phase_field_process_id_));

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    GlobalExecutor::executeSelectedMemberDereferenced(
        global_assembler_, &VectorMatrixAssembler::assembleWithJacobian,
        local_assemblers_, pv.getActiveElementIDs(), dof_tables, t, dt, x, xdot,
        dxdot_dx, dx_dx, process_id, M, K, b, Jac, coupled_solutions_);
}

template <int DisplacementDim>
void ThermoMechanicalPhaseFieldProcess<DisplacementDim>::
    preTimestepConcreteProcess(std::vector<GlobalVector*> const& x,
                               double const t,
                               double const dt,
                               const int process_id)
{
    DBUG("PreTimestep ThermoMechanicalPhaseFieldProcess.");

    if (process_id != mechanics_related_process_id_)
    {
        return;
    }

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    GlobalExecutor::executeSelectedMemberOnDereferenced(
        &ThermoMechanicalPhaseFieldLocalAssemblerInterface::preTimestep,
        local_assemblers_, pv.getActiveElementIDs(), getDOFTable(process_id),
        *x[process_id], t, dt);
}

template <int DisplacementDim>
void ThermoMechanicalPhaseFieldProcess<DisplacementDim>::
    postTimestepConcreteProcess(std::vector<GlobalVector*> const& x,
                                double const t,
                                double const dt,
                                int const process_id)
{
    DBUG("PostTimestep ThermoMechanicalPhaseFieldProcess.");

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    GlobalExecutor::executeSelectedMemberOnDereferenced(
        &ThermoMechanicalPhaseFieldLocalAssemblerInterface::postTimestep,
        local_assemblers_, pv.getActiveElementIDs(), getDOFTable(process_id),
        *x[process_id], t, dt);
}

template <int DisplacementDim>
void ThermoMechanicalPhaseFieldProcess<
    DisplacementDim>::postNonLinearSolverConcreteProcess(GlobalVector const& x,
                                                         const double t,
                                                         double const dt,
                                                         const int process_id)
{
    if (process_id != mechanics_related_process_id_)
    {
        return;
    }

    DBUG("PostNonLinearSolver ThermoMechanicalPhaseFieldProcess.");
    // Calculate strain, stress or other internal variables of mechanics.
    const bool use_monolithic_scheme = false;
    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    GlobalExecutor::executeSelectedMemberOnDereferenced(
        &LocalAssemblerInterface::postNonLinearSolver, local_assemblers_,
        pv.getActiveElementIDs(), getDOFTable(process_id), x, t, dt,
        use_monolithic_scheme);
}

template class ThermoMechanicalPhaseFieldProcess<2>;
template class ThermoMechanicalPhaseFieldProcess<3>;

}  // namespace ThermoMechanicalPhaseField
}  // namespace ProcessLib
