/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2023, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "PhaseFieldProcess.h"

#include <cassert>

#include "MeshLib/Utils/getOrCreateMeshProperty.h"
#include "NumLib/DOF/ComputeSparsityPattern.h"
#include "PhaseFieldFEM.h"
#include "ProcessLib/Process.h"
#include "ProcessLib/SmallDeformation/CreateLocalAssemblers.h"

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
      _process_data(std::move(process_data))
{
    if (use_monolithic_scheme)
    {
        OGS_FATAL(
            "Monolithic scheme is not implemented for the PhaseField process.");
    }

    _nodal_forces = MeshLib::getOrCreateMeshProperty<double>(
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
    // For the M process (deformation) in the staggered scheme.
    if (process_id == 0)
    {
        return *_local_to_global_index_map;
    }

    // For the equation of phasefield
    return *_local_to_global_index_map_single_component;
}

template <int DisplacementDim>
void PhaseFieldProcess<DisplacementDim>::constructDofTable()
{
    // For displacement equation.
    const int mechanics_process_id = 0;
    constructDofTableOfSpecifiedProcessStaggeredScheme(mechanics_process_id);

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

    // For phase field equation.
    _sparsity_pattern_with_single_component = NumLib::computeSparsityPattern(
        *_local_to_global_index_map_single_component, _mesh);
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
        NumLib::IntegrationOrder{integration_order}, mesh.isAxiallySymmetric(),
        _process_data);

    _secondary_variables.addSecondaryVariable(
        "sigma",
        makeExtrapolator(MathLib::KelvinVector::KelvinVectorType<
                             DisplacementDim>::RowsAtCompileTime,
                         getExtrapolator(), _local_assemblers,
                         &LocalAssemblerInterface::getIntPtSigma));

    _secondary_variables.addSecondaryVariable(
        "epsilon",
        makeExtrapolator(MathLib::KelvinVector::KelvinVectorType<
                             DisplacementDim>::RowsAtCompileTime,
                         getExtrapolator(), _local_assemblers,
                         &LocalAssemblerInterface::getIntPtEpsilon));

    _secondary_variables.addSecondaryVariable(
        "sigma_tensile",
        makeExtrapolator(MathLib::KelvinVector::KelvinVectorType<
                             DisplacementDim>::RowsAtCompileTime,
                         getExtrapolator(), _local_assemblers,
                         &LocalAssemblerInterface::getIntPtSigmaTensile));

    _secondary_variables.addSecondaryVariable(
        "sigma_compressive",
        makeExtrapolator(MathLib::KelvinVector::KelvinVectorType<
                             DisplacementDim>::RowsAtCompileTime,
                         getExtrapolator(), _local_assemblers,
                         &LocalAssemblerInterface::getIntPtSigmaCompressive));

    _secondary_variables.addSecondaryVariable(
        "eps_tensile",
        makeExtrapolator(MathLib::KelvinVector::KelvinVectorType<
                             DisplacementDim>::RowsAtCompileTime,
                         getExtrapolator(), _local_assemblers,
                         &LocalAssemblerInterface::getIntPtEpsilonTensile));

    // Initialize local assemblers after all variables have been set.
    GlobalExecutor::executeMemberOnDereferenced(
        &LocalAssemblerInterface::initialize, _local_assemblers,
        *_local_to_global_index_map);
}

template <int DisplacementDim>
void PhaseFieldProcess<DisplacementDim>::initializeBoundaryConditions(
    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const& media)
{
    // Staggered scheme:
    // for the equations of deformation.
    const int mechanical_process_id = 0;
    initializeProcessBoundaryConditionsAndSourceTerms(
        *_local_to_global_index_map, mechanical_process_id, media);
    // for the phase field
    const int phasefield_process_id = 1;
    initializeProcessBoundaryConditionsAndSourceTerms(
        *_local_to_global_index_map_single_component, phasefield_process_id,
        media);
}

template <int DisplacementDim>
void PhaseFieldProcess<DisplacementDim>::assembleConcreteProcess(
    const double t, double const dt, std::vector<GlobalVector*> const& x,
    std::vector<GlobalVector*> const& x_prev, int const process_id,
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
    dof_tables.emplace_back(*_local_to_global_index_map_single_component);
    dof_tables.emplace_back(*_local_to_global_index_map);

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    // Call global assembler for each local assembly item.
    GlobalExecutor::executeSelectedMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assemble, _local_assemblers,
        pv.getActiveElementIDs(), dof_tables, t, dt, x, x_prev, process_id, M,
        K, b);
}

template <int DisplacementDim>
void PhaseFieldProcess<DisplacementDim>::assembleWithJacobianConcreteProcess(
    const double t, double const dt, std::vector<GlobalVector*> const& x,
    std::vector<GlobalVector*> const& x_prev, int const process_id,
    GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b, GlobalMatrix& Jac)
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
    dof_tables.emplace_back(*_local_to_global_index_map);
    dof_tables.emplace_back(*_local_to_global_index_map_single_component);

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    // Call global assembler for each local assembly item.
    GlobalExecutor::executeSelectedMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assembleWithJacobian,
        _local_assemblers, pv.getActiveElementIDs(), dof_tables, t, dt, x,
        x_prev, process_id, M, K, b, Jac);

    if (process_id == 0)
    {
        b.copyValues(*_nodal_forces);
        std::transform(_nodal_forces->begin(), _nodal_forces->end(),
                       _nodal_forces->begin(), [](double val) { return -val; });
    }
}

template <int DisplacementDim>
void PhaseFieldProcess<DisplacementDim>::preTimestepConcreteProcess(
    std::vector<GlobalVector*> const& x, double const t, double const dt,
    const int process_id)
{
    DBUG("PreTimestep PhaseFieldProcess {:d}.", process_id);

    _process_data.injected_volume = t;

    _x_previous_timestep =
        MathLib::MatrixVectorTraits<GlobalVector>::newInstance(*x[process_id]);

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    GlobalExecutor::executeSelectedMemberOnDereferenced(
        &LocalAssemblerInterface::preTimestep, _local_assemblers,
        pv.getActiveElementIDs(), getDOFTable(process_id), *x[process_id], t,
        dt);
}

template <int DisplacementDim>
void PhaseFieldProcess<DisplacementDim>::postTimestepConcreteProcess(
    std::vector<GlobalVector*> const& x,
    std::vector<GlobalVector*> const& /*x_prev*/, const double t,
    const double /*delta_t*/, int const process_id)
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

        ProcessLib::ProcessVariable const& pv =
            getProcessVariables(process_id)[0];

        GlobalExecutor::executeSelectedMemberOnDereferenced(
            &LocalAssemblerInterface::computeEnergy, _local_assemblers,
            pv.getActiveElementIDs(), dof_tables, x, t,
            _process_data.elastic_energy, _process_data.surface_energy,
            _process_data.pressure_work);

#ifdef USE_PETSC
        double const elastic_energy = _process_data.elastic_energy;
        MPI_Allreduce(&elastic_energy, &_process_data.elastic_energy, 1,
                      MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
        double const surface_energy = _process_data.surface_energy;
        MPI_Allreduce(&surface_energy, &_process_data.surface_energy, 1,
                      MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
        double const pressure_work = _process_data.pressure_work;
        MPI_Allreduce(&pressure_work, &_process_data.pressure_work, 1,
                      MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
#endif

        INFO(
            "Elastic energy: {:g} Surface energy: {:g} Pressure work: {:g}  at "
            "time: {:g} ",
            _process_data.elastic_energy, _process_data.surface_energy,
            _process_data.pressure_work, t);
        if (_process_data.propagating_pressurized_crack)
        {
            INFO("Pressure: {:g} at time: {:g} ", _process_data.pressure, t);
        }
    }
}

template <int DisplacementDim>
void PhaseFieldProcess<DisplacementDim>::postNonLinearSolverConcreteProcess(
    std::vector<GlobalVector*> const& x,
    std::vector<GlobalVector*> const& /*x_prev*/, const double t,
    double const /*dt*/, const int process_id)
{
    _process_data.crack_volume = 0.0;

    if (isPhaseFieldProcess(process_id))
    {
        if (_process_data.propagating_pressurized_crack)
        {
            auto& u = *x[0];
            MathLib::LinAlg::scale(const_cast<GlobalVector&>(u),
                                   1 / _process_data.pressure);
        }
        return;
    }

    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_tables;

    dof_tables.emplace_back(*_local_to_global_index_map);
    dof_tables.emplace_back(*_local_to_global_index_map_single_component);

    DBUG("PostNonLinearSolver crack volume computation.");

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];
    GlobalExecutor::executeSelectedMemberOnDereferenced(
        &LocalAssemblerInterface::computeCrackIntegral, _local_assemblers,
        pv.getActiveElementIDs(), dof_tables, x, t, _process_data.crack_volume);

#ifdef USE_PETSC
    double const crack_volume = _process_data.crack_volume;
    MPI_Allreduce(&crack_volume, &_process_data.crack_volume, 1, MPI_DOUBLE,
                  MPI_SUM, PETSC_COMM_WORLD);
#endif

    INFO("Integral of crack: {:g}", _process_data.crack_volume);

    if (_process_data.propagating_pressurized_crack)
    {
        _process_data.pressure_old = _process_data.pressure;
        _process_data.pressure =
            _process_data.injected_volume / _process_data.crack_volume;
        _process_data.pressure_error =
            std::abs(_process_data.pressure_old - _process_data.pressure) /
            _process_data.pressure;
        INFO("Internal pressure: {:g} and Pressure error: {:.4e}",
             _process_data.pressure, _process_data.pressure_error);

        auto& u = *x[0];
        MathLib::LinAlg::scale(const_cast<GlobalVector&>(u),
                               _process_data.pressure);
    }
}

template <int DisplacementDim>
void PhaseFieldProcess<DisplacementDim>::updateConstraints(
    GlobalVector& lower, GlobalVector& upper, int const /*process_id*/)
{
    lower.setZero();
    MathLib::LinAlg::setLocalAccessibleVector(*_x_previous_timestep);
    MathLib::LinAlg::copy(*_x_previous_timestep, upper);

    GlobalIndexType const x_begin = _x_previous_timestep->getRangeBegin();
    GlobalIndexType const x_end = _x_previous_timestep->getRangeEnd();

    for (GlobalIndexType i = x_begin; i < x_end; i++)
    {
        if ((*_x_previous_timestep)[i] > _process_data.irreversible_threshold)
        {
            upper.set(i, 1.0);
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
