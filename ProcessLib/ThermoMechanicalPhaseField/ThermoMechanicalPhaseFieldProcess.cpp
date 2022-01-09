/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ThermoMechanicalPhaseFieldProcess.h"

#include <cassert>

#include "NumLib/DOF/ComputeSparsityPattern.h"
#include "ProcessLib/Process.h"
#include "ProcessLib/SmallDeformation/CreateLocalAssemblers.h"
#include "ThermoMechanicalPhaseFieldFEM.h"
#include "ThermoMechanicalPhaseFieldProcessData.h"

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
ThermoMechanicalPhaseFieldProcess<DisplacementDim>::getDOFTable(
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
    // For displacement equation.
    constructDofTableOfSpecifiedProcessStaggeredScheme(
        _mechanics_related_process_id);

    // TODO move the two data members somewhere else.
    // for extrapolation of secondary variables of stress or strain
    std::vector<MeshLib::MeshSubset> all_mesh_subsets_single_component{
        *_mesh_subset_all_nodes};
    _local_to_global_index_map_single_component =
        std::make_unique<NumLib::LocalToGlobalIndexMap>(
            std::move(all_mesh_subsets_single_component),
            // by location order is needed for output
            NumLib::ComponentOrder::BY_LOCATION);

    assert(_local_to_global_index_map_single_component);

    // For phase field equation or the heat conduction.
    _sparsity_pattern_with_single_component = NumLib::computeSparsityPattern(
        *_local_to_global_index_map_single_component, _mesh);
}

template <int DisplacementDim>
void ThermoMechanicalPhaseFieldProcess<DisplacementDim>::
    initializeConcreteProcess(NumLib::LocalToGlobalIndexMap const& dof_table,
                              MeshLib::Mesh const& mesh,
                              unsigned const integration_order)
{
    ProcessLib::SmallDeformation::createLocalAssemblers<
        DisplacementDim, ThermoMechanicalPhaseFieldLocalAssembler>(
        mesh.getElements(), dof_table, _local_assemblers,
        mesh.isAxiallySymmetric(), integration_order, _process_data,
        _mechanics_related_process_id, _phase_field_process_id,
        _heat_conduction_process_id);

    _secondary_variables.addSecondaryVariable(
        "sigma",
        makeExtrapolator(
            MathLib::KelvinVector::KelvinVectorType<
                DisplacementDim>::RowsAtCompileTime,
            getExtrapolator(), _local_assemblers,
            &ThermoMechanicalPhaseFieldLocalAssemblerInterface::getIntPtSigma));

    _secondary_variables.addSecondaryVariable(
        "epsilon",
        makeExtrapolator(MathLib::KelvinVector::KelvinVectorType<
                             DisplacementDim>::RowsAtCompileTime,
                         getExtrapolator(), _local_assemblers,
                         &ThermoMechanicalPhaseFieldLocalAssemblerInterface::
                             getIntPtEpsilon));

    _secondary_variables.addSecondaryVariable(
        "heat_flux",
        makeExtrapolator(mesh.getDimension(), getExtrapolator(),
                         _local_assemblers,
                         &ThermoMechanicalPhaseFieldLocalAssemblerInterface::
                             getIntPtHeatFlux));

    // Initialize local assemblers after all variables have been set.
    GlobalExecutor::executeMemberOnDereferenced(
        &LocalAssemblerInterface::initialize, _local_assemblers,
        *_local_to_global_index_map);
}

template <int DisplacementDim>
void ThermoMechanicalPhaseFieldProcess<
    DisplacementDim>::initializeBoundaryConditions()
{
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
void ThermoMechanicalPhaseFieldProcess<DisplacementDim>::
    assembleConcreteProcess(const double t, double const dt,
                            std::vector<GlobalVector*> const& x,
                            std::vector<GlobalVector*> const& xdot,
                            int const process_id, GlobalMatrix& M,
                            GlobalMatrix& K, GlobalVector& b)
{
    DBUG("Assemble the equations for ThermoMechanicalPhaseFieldProcess.");

    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_table = {std::ref(*_local_to_global_index_map)};
    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    // Call global assembler for each local assembly item.
    GlobalExecutor::executeSelectedMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assemble, _local_assemblers,
        pv.getActiveElementIDs(), dof_table, t, dt, x, xdot, process_id, M, K,
        b);
}

template <int DisplacementDim>
void ThermoMechanicalPhaseFieldProcess<DisplacementDim>::
    assembleWithJacobianConcreteProcess(const double t, double const dt,
                                        std::vector<GlobalVector*> const& x,
                                        std::vector<GlobalVector*> const& xdot,
                                        int const process_id, GlobalMatrix& M,
                                        GlobalMatrix& K, GlobalVector& b,
                                        GlobalMatrix& Jac)
{
    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_tables;
    // For the staggered scheme
    if (process_id == _mechanics_related_process_id)
    {
        DBUG(
            "Assemble the Jacobian equations of "
            "temperature-deformation in ThermoMechanicalPhaseFieldProcess for "
            "the staggered scheme.");
    }

    if (process_id == _phase_field_process_id)
    {
        DBUG(
            "Assemble the Jacobian equations ofphase field in "
            "ThermoMechanicalPhaseFieldProcess for the staggered scheme.");
    }
    else
    {
        DBUG(
            "Assemble the Jacobian equations of heat conduction in "
            "ThermoMechanicalPhaseFieldProcess for the staggered scheme.");
    }
    dof_tables.emplace_back(
        getDOFTableByProcessID(_heat_conduction_process_id));
    dof_tables.emplace_back(
        getDOFTableByProcessID(_mechanics_related_process_id));
    dof_tables.emplace_back(getDOFTableByProcessID(_phase_field_process_id));

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    GlobalExecutor::executeSelectedMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assembleWithJacobian,
        _local_assemblers, pv.getActiveElementIDs(), dof_tables, t, dt, x, xdot,
        process_id, M, K, b, Jac);
}

template <int DisplacementDim>
void ThermoMechanicalPhaseFieldProcess<DisplacementDim>::
    preTimestepConcreteProcess(std::vector<GlobalVector*> const& x,
                               double const t,
                               double const dt,
                               const int process_id)
{
    DBUG("PreTimestep ThermoMechanicalPhaseFieldProcess.");

    if (process_id != _mechanics_related_process_id)
    {
        return;
    }

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    GlobalExecutor::executeSelectedMemberOnDereferenced(
        &ThermoMechanicalPhaseFieldLocalAssemblerInterface::preTimestep,
        _local_assemblers, pv.getActiveElementIDs(), getDOFTable(process_id),
        *x[process_id], t, dt);
}

template <int DisplacementDim>
void ThermoMechanicalPhaseFieldProcess<DisplacementDim>::
    postTimestepConcreteProcess(std::vector<GlobalVector*> const& x,
                                double const t,
                                double const dt,
                                int const process_id)
{
    if (process_id != 0)
    {
        return;
    }

    DBUG("PostTimestep ThermoMechanicalPhaseFieldProcess.");
    std::vector<NumLib::LocalToGlobalIndexMap const*> dof_tables;
    auto const n_processes = x.size();
    dof_tables.reserve(n_processes);
    for (std::size_t process_id = 0; process_id < n_processes; ++process_id)
    {
        dof_tables.push_back(&getDOFTable(process_id));
    }

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];
    GlobalExecutor::executeSelectedMemberOnDereferenced(
        &ThermoMechanicalPhaseFieldLocalAssemblerInterface::postTimestep,
        _local_assemblers, pv.getActiveElementIDs(), dof_tables, x, t, dt);
}

template <int DisplacementDim>
void ThermoMechanicalPhaseFieldProcess<DisplacementDim>::
    postNonLinearSolverConcreteProcess(GlobalVector const& x,
                                       GlobalVector const& xdot, const double t,
                                       double const dt, const int process_id)
{
    if (process_id != _mechanics_related_process_id)
    {
        return;
    }

    DBUG("PostNonLinearSolver ThermoMechanicalPhaseFieldProcess.");
    // Calculate strain, stress or other internal variables of mechanics.
    const bool use_monolithic_scheme = false;
    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    GlobalExecutor::executeSelectedMemberOnDereferenced(
        &LocalAssemblerInterface::postNonLinearSolver, _local_assemblers,
        pv.getActiveElementIDs(), getDOFTable(process_id), x, xdot, t, dt,
        use_monolithic_scheme, process_id);
}

template class ThermoMechanicalPhaseFieldProcess<2>;
template class ThermoMechanicalPhaseFieldProcess<3>;

}  // namespace ThermoMechanicalPhaseField
}  // namespace ProcessLib
