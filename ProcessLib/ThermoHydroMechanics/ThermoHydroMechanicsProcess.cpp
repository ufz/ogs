/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ThermoHydroMechanicsProcess.h"

#include <cassert>

#include "MeshLib/Elements/Utils.h"
#include "NumLib/DOF/ComputeSparsityPattern.h"
#include "ProcessLib/Process.h"
#include "ProcessLib/ThermoHydroMechanics/CreateLocalAssemblers.h"

#include "ThermoHydroMechanicsFEM.h"
#include "ThermoHydroMechanicsProcessData.h"

namespace ProcessLib
{
namespace ThermoHydroMechanics
{
template <int DisplacementDim>
ThermoHydroMechanicsProcess<DisplacementDim>::ThermoHydroMechanicsProcess(
    std::string name,
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    unsigned const integration_order,
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
        process_variables,
    ThermoHydroMechanicsProcessData<DisplacementDim>&& process_data,
    SecondaryVariableCollection&& secondary_variables,
    bool const use_monolithic_scheme)
    : Process(std::move(name), mesh, std::move(jacobian_assembler), parameters,
              integration_order, std::move(process_variables),
              std::move(secondary_variables), use_monolithic_scheme),
      process_data_(std::move(process_data))
{
    nodal_forces_ = MeshLib::getOrCreateMeshProperty<double>(
        mesh, "NodalForces", MeshLib::MeshItemType::Node, DisplacementDim);

    hydraulic_flow_ = MeshLib::getOrCreateMeshProperty<double>(
        mesh, "HydraulicFlow", MeshLib::MeshItemType::Node, 1);
}

template <int DisplacementDim>
bool ThermoHydroMechanicsProcess<DisplacementDim>::isLinear() const
{
    return false;
}

template <int DisplacementDim>
MathLib::MatrixSpecifications
ThermoHydroMechanicsProcess<DisplacementDim>::getMatrixSpecifications(
    const int process_id) const
{
    // For the monolithic scheme or the M process (deformation) in the staggered
    // scheme.
    if (use_monolithic_scheme_ || process_id == 2)
    {
        auto const& l = *local_to_global_index_map_;
        return {l.dofSizeWithoutGhosts(), l.dofSizeWithoutGhosts(),
                &l.getGhostIndices(), &this->sparsity_pattern_};
    }

    // For staggered scheme and T or H process (pressure).
    auto const& l = *local_to_global_index_map_with_base_nodes_;
    return {l.dofSizeWithoutGhosts(), l.dofSizeWithoutGhosts(),
            &l.getGhostIndices(), &sparsity_pattern_with_linear_element_};
}

template <int DisplacementDim>
void ThermoHydroMechanicsProcess<DisplacementDim>::constructDofTable()
{
    // Create single component dof in every of the mesh's nodes.
    mesh_subset_all_nodes_ =
        std::make_unique<MeshLib::MeshSubset>(mesh_, mesh_.getNodes());
    // Create single component dof in the mesh's base nodes.
    base_nodes_ = MeshLib::getBaseNodes(mesh_.getElements());
    mesh_subset_base_nodes_ =
        std::make_unique<MeshLib::MeshSubset>(mesh_, base_nodes_);

    // TODO move the two data members somewhere else.
    // for extrapolation of secondary variables of stress or strain
    std::vector<MeshLib::MeshSubset> all_mesh_subsets_single_component{
        *mesh_subset_all_nodes_};
    local_to_global_index_map_single_component_ =
        std::make_unique<NumLib::LocalToGlobalIndexMap>(
            std::move(all_mesh_subsets_single_component),
            // by location order is needed for output
            NumLib::ComponentOrder::BY_LOCATION);

    if (use_monolithic_scheme_)
    {
        // For temperature, which is the first
        std::vector<MeshLib::MeshSubset> all_mesh_subsets{
            *mesh_subset_base_nodes_};

        // For pressure, which is the second
        all_mesh_subsets.push_back(*mesh_subset_base_nodes_);

        // For displacement.
        const int monolithic_process_id = 0;
        std::generate_n(std::back_inserter(all_mesh_subsets),
                        getProcessVariables(monolithic_process_id)[2]
                            .get()
                            .getNumberOfComponents(),
                        [&]() { return *mesh_subset_all_nodes_; });

        std::vector<int> const vec_n_components{1, 1, DisplacementDim};
        local_to_global_index_map_ =
            std::make_unique<NumLib::LocalToGlobalIndexMap>(
                std::move(all_mesh_subsets), vec_n_components,
                NumLib::ComponentOrder::BY_LOCATION);
        assert(local_to_global_index_map_);
    }
    else
    {
        // For displacement equation.
        const int process_id = 2;
        std::vector<MeshLib::MeshSubset> all_mesh_subsets;
        std::generate_n(
            std::back_inserter(all_mesh_subsets),
            getProcessVariables(process_id)[0].get().getNumberOfComponents(),
            [&]() { return *mesh_subset_all_nodes_; });

        std::vector<int> const vec_n_components{DisplacementDim};
        local_to_global_index_map_ =
            std::make_unique<NumLib::LocalToGlobalIndexMap>(
                std::move(all_mesh_subsets), vec_n_components,
                NumLib::ComponentOrder::BY_LOCATION);

        // For pressure equation or temperature equation.
        // Collect the mesh subsets with base nodes in a vector.
        std::vector<MeshLib::MeshSubset> all_mesh_subsets_base_nodes{
            *mesh_subset_base_nodes_};
        local_to_global_index_map_with_base_nodes_ =
            std::make_unique<NumLib::LocalToGlobalIndexMap>(
                std::move(all_mesh_subsets_base_nodes),
                // by location order is needed for output
                NumLib::ComponentOrder::BY_LOCATION);

        sparsity_pattern_with_linear_element_ = NumLib::computeSparsityPattern(
            *local_to_global_index_map_with_base_nodes_, mesh_);

        assert(local_to_global_index_map_);
        assert(local_to_global_index_map_with_base_nodes_);
    }
}

template <int DisplacementDim>
void ThermoHydroMechanicsProcess<DisplacementDim>::initializeConcreteProcess(
    NumLib::LocalToGlobalIndexMap const& dof_table,
    MeshLib::Mesh const& mesh,
    unsigned const integration_order)
{
    const int mechanical_process_id = use_monolithic_scheme_ ? 0 : 2;
    const int deformation_variable_id = use_monolithic_scheme_ ? 2 : 0;
    ProcessLib::ThermoHydroMechanics::createLocalAssemblers<
        DisplacementDim, ThermoHydroMechanicsLocalAssembler>(
        mesh.getDimension(), mesh.getElements(), dof_table,
        // use displacement process variable to set shape function order
        getProcessVariables(mechanical_process_id)[deformation_variable_id]
            .get()
            .getShapeFunctionOrder(),
        local_assemblers_, mesh.isAxiallySymmetric(), integration_order,
        process_data_);

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

    secondary_variables_.addSecondaryVariable(
        "velocity",
        makeExtrapolator(mesh.getDimension(), getExtrapolator(),
                         local_assemblers_,
                         &LocalAssemblerInterface::getIntPtDarcyVelocity));

    process_data_.pressure_interpolated =
        MeshLib::getOrCreateMeshProperty<double>(
            const_cast<MeshLib::Mesh&>(mesh), "pressure_interpolated",
            MeshLib::MeshItemType::Node, 1);

    process_data_.temperature_interpolated =
        MeshLib::getOrCreateMeshProperty<double>(
            const_cast<MeshLib::Mesh&>(mesh), "temperature_interpolated",
            MeshLib::MeshItemType::Node, 1);

    // Initialize local assemblers after all variables have been set.
    GlobalExecutor::executeMemberOnDereferenced(
        &LocalAssemblerInterface::initialize, local_assemblers_,
        *local_to_global_index_map_);
}

template <int DisplacementDim>
void ThermoHydroMechanicsProcess<
    DisplacementDim>::initializeBoundaryConditions()
{
    if (use_monolithic_scheme_)
    {
        const int process_id_of_thermohydromechancs = 0;
        initializeProcessBoundaryConditionsAndSourceTerms(
            *local_to_global_index_map_, process_id_of_thermohydromechancs);
        return;
    }

    // Staggered scheme:
    // for the equations of heat transport
    const int thermal_process_id = 0;
    initializeProcessBoundaryConditionsAndSourceTerms(
        *local_to_global_index_map_with_base_nodes_, thermal_process_id);

    // for the equations of mass balance
    const int hydraulic_process_id = 1;
    initializeProcessBoundaryConditionsAndSourceTerms(
        *local_to_global_index_map_with_base_nodes_, hydraulic_process_id);

    // for the equations of deformation.
    const int mechanical_process_id = 2;
    initializeProcessBoundaryConditionsAndSourceTerms(
        *local_to_global_index_map_, mechanical_process_id);
}

template <int DisplacementDim>
void ThermoHydroMechanicsProcess<DisplacementDim>::assembleConcreteProcess(
    const double t, double const dt, std::vector<GlobalVector*> const& x,
    std::vector<GlobalVector*> const& xdot, int const process_id,
    GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b)
{
    DBUG("Assemble the equations for ThermoHydroMechanics");

    // Note: This assembly function is for the Picard nonlinear solver. Since
    // only the Newton-Raphson method is employed to simulate coupled HM
    // processes in this class, this function is actually not used so far.

    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_table = {std::ref(*local_to_global_index_map_)};
    // Call global assembler for each local assembly item.
    GlobalExecutor::executeMemberDereferenced(
        global_assembler_, &VectorMatrixAssembler::assemble, local_assemblers_,
        dof_table, t, dt, x, xdot, process_id, M, K, b, coupled_solutions_);
}

template <int DisplacementDim>
void ThermoHydroMechanicsProcess<DisplacementDim>::
    assembleWithJacobianConcreteProcess(
        const double t, double const dt, std::vector<GlobalVector*> const& x,
        GlobalVector const& xdot, const double dxdot_dx, const double dx_dx,
        int const process_id, GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b,
        GlobalMatrix& Jac)
{
    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_tables;
    // For the monolithic scheme
    if (use_monolithic_scheme_)
    {
        DBUG(
            "Assemble the Jacobian of ThermoHydroMechanics for the monolithic"
            " scheme.");
        dof_tables.emplace_back(*local_to_global_index_map_);
    }
    else
    {
        // For the staggered scheme
        if (process_id == 0)
        {
            DBUG(
                "Assemble the Jacobian equations of heat transport process in "
                "ThermoHydroMechanics for the staggered scheme.");
        }
        else if (process_id == 1)
        {
            DBUG(
                "Assemble the Jacobian equations of liquid fluid process in "
                "ThermoHydroMechanics for the staggered scheme.");
        }
        else
        {
            DBUG(
                "Assemble the Jacobian equations of mechanical process in "
                "ThermoHydroMechanics for the staggered scheme.");
        }
        dof_tables.emplace_back(*local_to_global_index_map_with_base_nodes_);
        dof_tables.emplace_back(*local_to_global_index_map_with_base_nodes_);
        dof_tables.emplace_back(*local_to_global_index_map_);
    }

    GlobalExecutor::executeMemberDereferenced(
        global_assembler_, &VectorMatrixAssembler::assembleWithJacobian,
        local_assemblers_, dof_tables, t, dt, x, xdot, dxdot_dx, dx_dx,
        process_id, M, K, b, Jac, coupled_solutions_);

    auto copyRhs = [&](int const variable_id, auto& output_vector) {
        if (use_monolithic_scheme_)
        {
            transformVariableFromGlobalVector(b, variable_id, dof_tables[0],
                                              output_vector,
                                              std::negate<double>());
        }
        else
        {
            transformVariableFromGlobalVector(b, 0, dof_tables[process_id],
                                              output_vector,
                                              std::negate<double>());
        }
    };
    if (use_monolithic_scheme_ || process_id == 1)
    {
        copyRhs(0, *hydraulic_flow_);
    }
    if (use_monolithic_scheme_ || process_id == 2)
    {
        copyRhs(1, *nodal_forces_);
    }
}

template <int DisplacementDim>
void ThermoHydroMechanicsProcess<DisplacementDim>::preTimestepConcreteProcess(
    std::vector<GlobalVector*> const& x, double const t, double const dt,
    const int process_id)
{
    DBUG("PreTimestep ThermoHydroMechanicsProcess.");

    if (hasMechanicalProcess(process_id))
    {
        GlobalExecutor::executeMemberOnDereferenced(
            &LocalAssemblerInterface::preTimestep, local_assemblers_,
            *local_to_global_index_map_, *x[process_id], t, dt);
    }
}

template <int DisplacementDim>
void ThermoHydroMechanicsProcess<DisplacementDim>::postTimestepConcreteProcess(
    std::vector<GlobalVector*> const& x, double const t, double const dt,
    const int process_id)
{
    DBUG("PostTimestep ThermoHydroMechanicsProcess.");
    GlobalExecutor::executeMemberOnDereferenced(
        &LocalAssemblerInterface::postTimestep, local_assemblers_,
        getDOFTable(process_id), *x[process_id], t, dt);
}

template <int DisplacementDim>
void ThermoHydroMechanicsProcess<
    DisplacementDim>::postNonLinearSolverConcreteProcess(GlobalVector const& x,
                                                         const double t,
                                                         double const dt,
                                                         const int process_id)
{
    if (!hasMechanicalProcess(process_id))
    {
        return;
    }

    DBUG("PostNonLinearSolver ThermoHydroMechanicsProcess.");
    // Calculate strain, stress or other internal variables of mechanics.
    GlobalExecutor::executeMemberOnDereferenced(
        &LocalAssemblerInterface::postNonLinearSolver, local_assemblers_,
        getDOFTable(process_id), x, t, dt, use_monolithic_scheme_);
}

template <int DisplacementDim>
void ThermoHydroMechanicsProcess<DisplacementDim>::
    computeSecondaryVariableConcrete(double const t, double const dt,
                                     GlobalVector const& x,
                                     GlobalVector const& x_dot,
                                     const int process_id)
{
    DBUG("Compute the secondary variables for ThermoHydroMechanicsProcess.");
    GlobalExecutor::executeMemberOnDereferenced(
        &LocalAssemblerInterface::computeSecondaryVariable, local_assemblers_,
        getDOFTable(process_id), t, dt, x, x_dot, coupled_solutions_);
}

template <int DisplacementDim>
std::tuple<NumLib::LocalToGlobalIndexMap*, bool> ThermoHydroMechanicsProcess<
    DisplacementDim>::getDOFTableForExtrapolatorData() const
{
    const bool manage_storage = false;
    return std::make_tuple(local_to_global_index_map_single_component_.get(),
                           manage_storage);
}

template <int DisplacementDim>
NumLib::LocalToGlobalIndexMap const&
ThermoHydroMechanicsProcess<DisplacementDim>::getDOFTable(
    const int process_id) const
{
    if (hasMechanicalProcess(process_id))
    {
        return *local_to_global_index_map_;
    }

    // For the equation of pressure
    return *local_to_global_index_map_with_base_nodes_;
}

template class ThermoHydroMechanicsProcess<2>;
template class ThermoHydroMechanicsProcess<3>;

}  // namespace ThermoHydroMechanics
}  // namespace ProcessLib
