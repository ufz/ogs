// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "HMPhaseFieldProcess.h"

#include <cassert>

#include "HMPhaseFieldFEM.h"
#include "MeshLib/Utils/getOrCreateMeshProperty.h"
#include "NumLib/DOF/ComputeSparsityPattern.h"
#include "ProcessLib/Process.h"
#include "ProcessLib/SmallDeformation/CreateLocalAssemblers.h"

namespace ProcessLib
{
namespace HMPhaseField
{
template <int DisplacementDim>
HMPhaseFieldProcess<DisplacementDim>::HMPhaseFieldProcess(
    std::string name,
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    unsigned const integration_order,
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
        process_variables,
    HMPhaseFieldProcessData<DisplacementDim>&& process_data,
    SecondaryVariableCollection&& secondary_variables,
    bool const use_monolithic_scheme)
    : Process(std::move(name), mesh, std::move(jacobian_assembler), parameters,
              integration_order, std::move(process_variables),
              std::move(secondary_variables), use_monolithic_scheme),
      _process_data(std::move(process_data))
{
    // For numerical Jacobian assembler
    if (this->_jacobian_assembler->isPerturbationEnabled())
    {
        OGS_FATAL(
            "The numericial Jacobian assembler is not supported for the "
            "HMPhaseField process.");
    }

    if (use_monolithic_scheme)
    {
        OGS_FATAL(
            "Monolithic scheme is not implemented for the HMPhaseField "
            "process.");
    }

    _nodal_forces = MeshLib::getOrCreateMeshProperty<double>(
        mesh, "NodalForces", MeshLib::MeshItemType::Node, DisplacementDim);
}

template <int DisplacementDim>
bool HMPhaseFieldProcess<DisplacementDim>::isLinear() const
{
    return false;
}

template <int DisplacementDim>
MathLib::MatrixSpecifications
HMPhaseFieldProcess<DisplacementDim>::getMatrixSpecifications(
    const int process_id) const
{
    // For the M process (deformation) in the staggered scheme.
    if (process_id == _process_data._mechanics_related_process_id)
    {
        auto const& l = *_local_to_global_index_map;
        return {l.dofSizeWithoutGhosts(), l.dofSizeWithoutGhosts(),
                &l.getGhostIndices(), &this->_sparsity_pattern};
    }

    // For the H (hydro) or PF (phasefield) process in the staggered scheme.
    auto const& l = *_local_to_global_index_map_single_component;
    return {l.dofSizeWithoutGhosts(), l.dofSizeWithoutGhosts(),
            &l.getGhostIndices(), &_sparsity_pattern_with_single_component};
}

template <int DisplacementDim>
NumLib::LocalToGlobalIndexMap const&
HMPhaseFieldProcess<DisplacementDim>::getDOFTable(const int process_id) const
{
    // For the M process (deformation) in the staggered scheme.
    if (process_id == _process_data._mechanics_related_process_id)
    {
        return *_local_to_global_index_map;
    }

    // For the H (hydro) or PF (phasefield) process in the staggered scheme.
    return *_local_to_global_index_map_single_component;
}

template <int DisplacementDim>
void HMPhaseFieldProcess<DisplacementDim>::constructDofTable()
{
    // For displacement equation.
    constructDofTableOfSpecifiedProcessStaggeredScheme(
        _process_data._mechanics_related_process_id);

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

    // For the H (hydro) or PF (phasefield) process in the staggered scheme.
    _sparsity_pattern_with_single_component = NumLib::computeSparsityPattern(
        *_local_to_global_index_map_single_component, _mesh);
}

template <int DisplacementDim>
void HMPhaseFieldProcess<DisplacementDim>::initializeConcreteProcess(
    NumLib::LocalToGlobalIndexMap const& dof_table,
    MeshLib::Mesh const& mesh,
    unsigned const integration_order)
{
    ProcessLib::SmallDeformation::createLocalAssemblers<
        DisplacementDim, HMPhaseFieldLocalAssembler>(
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
        "width_node",
        makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                         &LocalAssemblerInterface::getIntPtWidth));

    _process_data.ele_d = MeshLib::getOrCreateMeshProperty<double>(
        const_cast<MeshLib::Mesh&>(mesh), "damage", MeshLib::MeshItemType::Cell,
        1);

    _process_data.width = MeshLib::getOrCreateMeshProperty<double>(
        const_cast<MeshLib::Mesh&>(mesh), "width", MeshLib::MeshItemType::Cell,
        1);

    // Initialize local assemblers after all variables have been set.
    GlobalExecutor::executeMemberOnDereferenced(
        &LocalAssemblerInterface::initialize, _local_assemblers,
        *_local_to_global_index_map);
}

template <int DisplacementDim>
void HMPhaseFieldProcess<DisplacementDim>::initializeBoundaryConditions(
    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const& media)
{
    // Staggered scheme:
    // for the phase field
    initializeProcessBoundaryConditionsAndSourceTerms(
        *_local_to_global_index_map_single_component,
        _process_data._phasefield_process_id, media);
    // for the pressure
    initializeProcessBoundaryConditionsAndSourceTerms(
        *_local_to_global_index_map_single_component,
        _process_data._hydro_process_id, media);
    // for the deformation.
    initializeProcessBoundaryConditionsAndSourceTerms(
        *_local_to_global_index_map,
        _process_data._mechanics_related_process_id, media);
}

template <int DisplacementDim>
void HMPhaseFieldProcess<DisplacementDim>::assembleConcreteProcess(
    const double /*t*/, double const /*dt*/,
    std::vector<GlobalVector*> const& /*x*/,
    std::vector<GlobalVector*> const& /*x_prev*/, int const /*process_id*/,
    GlobalMatrix& /*M*/, GlobalMatrix& /*K*/, GlobalVector& /*b*/)
{
    OGS_FATAL(
        "HMPhaseFieldLocalAssembler: assembly for the Picard non-linear solver "
        "is not implemented.");
}

template <int DisplacementDim>
void HMPhaseFieldProcess<DisplacementDim>::assembleWithJacobianConcreteProcess(
    const double t, double const dt, std::vector<GlobalVector*> const& x,
    std::vector<GlobalVector*> const& x_prev, int const process_id,
    GlobalVector& b, GlobalMatrix& Jac)
{
    std::vector<NumLib::LocalToGlobalIndexMap const*> dof_tables;

    // For the staggered scheme
    if (process_id == _process_data._phasefield_process_id)
    {
        DBUG(
            "Assemble the Jacobian equations of phase field in "
            "HMPhaseFieldProcess for the staggered scheme.");
    }
    else if (process_id == _process_data._hydro_process_id)
    {
        DBUG(
            "Assemble the Jacobian equations of pressure in "
            "HMPhaseFieldProcess for the staggered scheme.");
    }
    else
    {
        DBUG(
            "Assemble the Jacobian equations of deformation in "
            "HMPhaseFieldProcess for the staggered scheme.");
    }

    dof_tables.emplace_back(_local_to_global_index_map_single_component.get());
    dof_tables.emplace_back(_local_to_global_index_map_single_component.get());
    dof_tables.emplace_back(_local_to_global_index_map.get());

    // Call global assembler for each local assembly item.
    GlobalExecutor::executeSelectedMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assembleWithJacobian,
        _local_assemblers, getActiveElementIDs(), dof_tables, t, dt, x, x_prev,
        process_id, &b, &Jac);

    if (process_id == _process_data._mechanics_related_process_id)
    {
        b.copyValues(*_nodal_forces);
        std::transform(_nodal_forces->begin(), _nodal_forces->end(),
                       _nodal_forces->begin(), [](double val) { return -val; });
    }
}

template <int DisplacementDim>
void HMPhaseFieldProcess<DisplacementDim>::preTimestepConcreteProcess(
    std::vector<GlobalVector*> const& x, double const t, double const dt,
    const int process_id)
{
    DBUG("PreTimestep HMPhaseFieldProcess {}.", process_id);

    if (process_id == _process_data._phasefield_process_id)
    {
        DBUG("Store the value of phase field at previous time step.");
        _x_previous_timestep =
            MathLib::MatrixVectorTraits<GlobalVector>::newInstance(
                *x[process_id]);
    }

    GlobalExecutor::executeSelectedMemberOnDereferenced(
        &LocalAssemblerInterface::preTimestep, _local_assemblers,
        getActiveElementIDs(), getDOFTable(process_id), *x[process_id], t, dt);
}

template <int DisplacementDim>
void HMPhaseFieldProcess<DisplacementDim>::postTimestepConcreteProcess(
    std::vector<GlobalVector*> const& x,
    std::vector<GlobalVector*> const& x_prev, const double t, const double dt,
    int const process_id)
{
    if (process_id == _process_data._phasefield_process_id)
    {
        DBUG("PostTimestep HMPhaseFieldProcess.");

        _process_data.elastic_energy = 0.0;
        _process_data.surface_energy = 0.0;
        _process_data.pressure_work = 0.0;

        std::vector<NumLib::LocalToGlobalIndexMap const*> dof_tables;

        dof_tables.emplace_back(
            _local_to_global_index_map_single_component.get());
        dof_tables.emplace_back(
            _local_to_global_index_map_single_component.get());
        dof_tables.emplace_back(_local_to_global_index_map.get());

        GlobalExecutor::executeSelectedMemberOnDereferenced(
            &LocalAssemblerInterface::computeEnergy, _local_assemblers,
            getActiveElementIDs(), dof_tables, x, t,
            _process_data.elastic_energy, _process_data.surface_energy,
            _process_data.pressure_work);

        showEnergyAndWork(t, _process_data.elastic_energy,
                          _process_data.surface_energy,
                          _process_data.pressure_work);
    }

    GlobalExecutor::executeSelectedMemberOnDereferenced(
        &LocalAssemblerInterface::postTimestep, _local_assemblers,
        getActiveElementIDs(), getDOFTables(x.size()), x, x_prev, t, dt,
        process_id);
}

template <int DisplacementDim>
void HMPhaseFieldProcess<DisplacementDim>::postNonLinearSolverConcreteProcess(
    std::vector<GlobalVector*> const& x,
    std::vector<GlobalVector*> const& x_prev, const double t, double const dt,
    const int process_id)
{
    std::vector<NumLib::LocalToGlobalIndexMap const*> dof_tables;

    dof_tables.emplace_back(_local_to_global_index_map_single_component.get());
    dof_tables.emplace_back(_local_to_global_index_map_single_component.get());
    dof_tables.emplace_back(_local_to_global_index_map.get());

    if (process_id == _process_data._phasefield_process_id)
    {
        INFO("Update fracture width and porous properties");
        GlobalExecutor::executeMemberOnDereferenced(
            &LocalAssemblerInterface::approximateFractureWidth,
            _local_assemblers, dof_tables, x, t, dt);
    }

    if (process_id == _process_data._hydro_process_id)
    {
        INFO("PostNonLinearSolver for Hydro Process");

        GlobalExecutor::executeSelectedMemberOnDereferenced(
            &LocalAssemblerInterface::postNonLinearSolver, _local_assemblers,
            getActiveElementIDs(), getDOFTables(x.size()), x, x_prev, t, dt,
            process_id);
    }
}

template <int DisplacementDim>
void HMPhaseFieldProcess<DisplacementDim>::updateConstraints(
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

template class HMPhaseFieldProcess<2>;
template class HMPhaseFieldProcess<3>;

}  // namespace HMPhaseField
}  // namespace ProcessLib
