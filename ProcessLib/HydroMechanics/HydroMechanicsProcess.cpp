/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "HydroMechanicsProcess.h"

#include <cassert>

#include "HydroMechanicsFEM.h"
#include "HydroMechanicsProcessData.h"
#include "MeshLib/Elements/Utils.h"
#include "NumLib/DOF/ComputeSparsityPattern.h"
#include "ProcessLib/Deformation/SolidMaterialInternalToSecondaryVariables.h"
#include "ProcessLib/HydroMechanics/CreateLocalAssemblers.h"
#include "ProcessLib/Process.h"
namespace ProcessLib
{
namespace HydroMechanics
{
template <int DisplacementDim>
HydroMechanicsProcess<DisplacementDim>::HydroMechanicsProcess(
    std::string name,
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    unsigned const integration_order,
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
        process_variables,
    HydroMechanicsProcessData<DisplacementDim>&& process_data,
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

    integration_point_writer_.emplace_back(
        std::make_unique<IntegrationPointWriter>(
            "sigma_ip",
            static_cast<int>(mesh.getDimension() == 2 ? 4 : 6) /*n components*/,
            2 /*integration order*/, [this]() {
                // Result containing integration point data for each local
                // assembler.
                std::vector<std::vector<double>> result;
                result.resize(local_assemblers_.size());

                for (std::size_t i = 0; i < local_assemblers_.size(); ++i)
                {
                    auto const& local_asm = *local_assemblers_[i];

                    result[i] = local_asm.getSigma();
                }

                return result;
            }));

    integration_point_writer_.emplace_back(
        std::make_unique<IntegrationPointWriter>(
            "epsilon_ip",
            static_cast<int>(mesh.getDimension() == 2 ? 4 : 6) /*n components*/,
            2 /*integration order*/, [this]() {
                // Result containing integration point data for each local
                // assembler.
                std::vector<std::vector<double>> result;
                result.resize(local_assemblers_.size());

                for (std::size_t i = 0; i < local_assemblers_.size(); ++i)
                {
                    auto const& local_asm = *local_assemblers_[i];

                    result[i] = local_asm.getEpsilon();
                }

                return result;
            }));
}

template <int DisplacementDim>
bool HydroMechanicsProcess<DisplacementDim>::isLinear() const
{
    return false;
}

template <int DisplacementDim>
MathLib::MatrixSpecifications
HydroMechanicsProcess<DisplacementDim>::getMatrixSpecifications(
    const int process_id) const
{
    // For the monolithic scheme or the M process (deformation) in the staggered
    // scheme.
    if (use_monolithic_scheme_ || process_id == 1)
    {
        auto const& l = *local_to_global_index_map_;
        return {l.dofSizeWithoutGhosts(), l.dofSizeWithoutGhosts(),
                &l.getGhostIndices(), &this->sparsity_pattern_};
    }

    // For staggered scheme and H process (pressure).
    auto const& l = *local_to_global_index_map_with_base_nodes_;
    return {l.dofSizeWithoutGhosts(), l.dofSizeWithoutGhosts(),
            &l.getGhostIndices(), &sparsity_pattern_with_linear_element_};
}

template <int DisplacementDim>
void HydroMechanicsProcess<DisplacementDim>::constructDofTable()
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
        // For pressure, which is the first
        std::vector<MeshLib::MeshSubset> all_mesh_subsets{
            *mesh_subset_base_nodes_};

        // For displacement.
        const int monolithic_process_id = 0;
        std::generate_n(std::back_inserter(all_mesh_subsets),
                        getProcessVariables(monolithic_process_id)[1]
                            .get()
                            .getNumberOfComponents(),
                        [&]() { return *mesh_subset_all_nodes_; });

        std::vector<int> const vec_n_components{1, DisplacementDim};
        local_to_global_index_map_ =
            std::make_unique<NumLib::LocalToGlobalIndexMap>(
                std::move(all_mesh_subsets), vec_n_components,
                NumLib::ComponentOrder::BY_LOCATION);
        assert(local_to_global_index_map_);
    }
    else
    {
        // For displacement equation.
        const int process_id = 1;
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

        // For pressure equation.
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
void HydroMechanicsProcess<DisplacementDim>::initializeConcreteProcess(
    NumLib::LocalToGlobalIndexMap const& dof_table,
    MeshLib::Mesh const& mesh,
    unsigned const integration_order)
{
    const int mechanical_process_id = use_monolithic_scheme_ ? 0 : 1;
    const int deformation_variable_id = use_monolithic_scheme_ ? 1 : 0;
    ProcessLib::HydroMechanics::createLocalAssemblers<
        DisplacementDim, HydroMechanicsLocalAssembler>(
        mesh.getDimension(), mesh.getElements(), dof_table,
        // use displacement process variable to set shape function order
        getProcessVariables(mechanical_process_id)[deformation_variable_id]
            .get()
            .getShapeFunctionOrder(),
        local_assemblers_, mesh.isAxiallySymmetric(), integration_order,
        process_data_);

    auto add_secondary_variable = [&](std::string const& name,
                                      int const num_components,
                                      auto get_ip_values_function) {
        secondary_variables_.addSecondaryVariable(
            name,
            makeExtrapolator(num_components, getExtrapolator(),
                             local_assemblers_,
                             std::move(get_ip_values_function)));
    };

    add_secondary_variable(
        "sigma_xx", 1, &LocalAssemblerIF::getIntPtSigmaXX);

    add_secondary_variable(
        "sigma_yy", 1, &LocalAssemblerIF::getIntPtSigmaYY);

    add_secondary_variable(
        "sigma_zz", 1, &LocalAssemblerIF::getIntPtSigmaZZ);

    add_secondary_variable(
        "sigma_xy", 1, &LocalAssemblerIF::getIntPtSigmaXY);

    if (DisplacementDim == 3)
    {
        add_secondary_variable(
            "sigma_xz", 1, &LocalAssemblerIF::getIntPtSigmaXZ);

        add_secondary_variable(
            "sigma_yz", 1, &LocalAssemblerIF::getIntPtSigmaYZ);
    }

    add_secondary_variable(
        "epsilon_xx", 1, &LocalAssemblerIF::getIntPtEpsilonXX);

    add_secondary_variable(
        "epsilon_yy", 1, &LocalAssemblerIF::getIntPtEpsilonYY);

    add_secondary_variable(
        "epsilon_zz", 1, &LocalAssemblerIF::getIntPtEpsilonZZ);

    add_secondary_variable(
        "epsilon_xy", 1, &LocalAssemblerIF::getIntPtEpsilonXY);

    if (DisplacementDim == 3)
    {
        add_secondary_variable(
            "epsilon_xz", 1, &LocalAssemblerIF::getIntPtEpsilonXZ);

        add_secondary_variable(
            "epsilon_yz", 1, &LocalAssemblerIF::getIntPtEpsilonYZ);
    }

    add_secondary_variable("velocity",
                           DisplacementDim,
                           &LocalAssemblerIF::getIntPtDarcyVelocity);

    //
    // enable output of internal variables defined by material models
    //
    ProcessLib::Deformation::solidMaterialInternalToSecondaryVariables<
        LocalAssemblerIF>(process_data_.solid_materials,
                                 add_secondary_variable);

    process_data_.pressure_interpolated =
        MeshLib::getOrCreateMeshProperty<double>(
            const_cast<MeshLib::Mesh&>(mesh), "pressure_interpolated",
            MeshLib::MeshItemType::Node, 1);

    process_data_.principal_stress_vector[0] =
        MeshLib::getOrCreateMeshProperty<double>(
            const_cast<MeshLib::Mesh&>(mesh), "principal_stress_vector_1",
            MeshLib::MeshItemType::Cell, 3);

    process_data_.principal_stress_vector[1] =
        MeshLib::getOrCreateMeshProperty<double>(
            const_cast<MeshLib::Mesh&>(mesh), "principal_stress_vector_2",
            MeshLib::MeshItemType::Cell, 3);

    process_data_.principal_stress_vector[2] =
        MeshLib::getOrCreateMeshProperty<double>(
            const_cast<MeshLib::Mesh&>(mesh), "principal_stress_vector_3",
            MeshLib::MeshItemType::Cell, 3);

    process_data_.principal_stress_values =
        MeshLib::getOrCreateMeshProperty<double>(
            const_cast<MeshLib::Mesh&>(mesh), "principal_stress_values",
            MeshLib::MeshItemType::Cell, 3);

    process_data_.permeability = MeshLib::getOrCreateMeshProperty<double>(
        const_cast<MeshLib::Mesh&>(mesh), "permeability",
        MeshLib::MeshItemType::Cell,
        MathLib::KelvinVector::KelvinVectorDimensions<DisplacementDim>::value);

    // Set initial conditions for integration point data.
    for (auto const& ip_writer : integration_point_writer_)
    {
        // Find the mesh property with integration point writer's name.
        auto const& name = ip_writer->name();
        if (!mesh.getProperties().existsPropertyVector<double>(name))
        {
            continue;
        }
        auto const& mesh_property =
            *mesh.getProperties().template getPropertyVector<double>(name);

        // The mesh property must be defined on integration points.
        if (mesh_property.getMeshItemType() !=
            MeshLib::MeshItemType::IntegrationPoint)
        {
            continue;
        }

        auto const ip_meta_data = getIntegrationPointMetaData(mesh, name);

        // Check the number of components.
        if (ip_meta_data.n_components != mesh_property.getNumberOfComponents())
        {
            OGS_FATAL(
                "Different number of components in meta data ({:d}) than in "
                "the integration point field data for '{:s}': {:d}.",
                ip_meta_data.n_components, name,
                mesh_property.getNumberOfComponents());
        }

        // Now we have a properly named vtk's field data array and the
        // corresponding meta data.
        std::size_t position = 0;
        for (auto& local_asm : local_assemblers_)
        {
            std::size_t const integration_points_read =
                local_asm->setIPDataInitialConditions(
                    name, &mesh_property[position],
                    ip_meta_data.integration_order);
            if (integration_points_read == 0)
            {
                OGS_FATAL(
                    "No integration points read in the integration point "
                    "initial conditions set function.");
            }
            position += integration_points_read * ip_meta_data.n_components;
        }
    }

    // Initialize local assemblers after all variables have been set.
    GlobalExecutor::executeMemberOnDereferenced(
        &LocalAssemblerIF::initialize, local_assemblers_,
        *local_to_global_index_map_);
}

template <int DisplacementDim>
void HydroMechanicsProcess<DisplacementDim>::initializeBoundaryConditions()
{
    if (use_monolithic_scheme_)
    {
        const int process_id_of_hydromechancs = 0;
        initializeProcessBoundaryConditionsAndSourceTerms(
            *local_to_global_index_map_, process_id_of_hydromechancs);
        return;
    }

    // Staggered scheme:
    // for the equations of pressure
    const int hydraulic_process_id = 0;
    initializeProcessBoundaryConditionsAndSourceTerms(
        *local_to_global_index_map_with_base_nodes_, hydraulic_process_id);

    // for the equations of deformation.
    const int mechanical_process_id = 1;
    initializeProcessBoundaryConditionsAndSourceTerms(
        *local_to_global_index_map_, mechanical_process_id);
}

template <int DisplacementDim>
void HydroMechanicsProcess<DisplacementDim>::assembleConcreteProcess(
    const double t, double const dt, std::vector<GlobalVector*> const& x,
    std::vector<GlobalVector*> const& xdot, int const process_id,
    GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b)
{
    DBUG("Assemble the equations for HydroMechanics");

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
void HydroMechanicsProcess<DisplacementDim>::
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
            "Assemble the Jacobian of HydroMechanics for the monolithic"
            " scheme.");
        dof_tables.emplace_back(*local_to_global_index_map_);
    }
    else
    {
        // For the staggered scheme
        if (process_id == 0)
        {
            DBUG(
                "Assemble the Jacobian equations of liquid fluid process in "
                "HydroMechanics for the staggered scheme.");
        }
        else
        {
            DBUG(
                "Assemble the Jacobian equations of mechanical process in "
                "HydroMechanics for the staggered scheme.");
        }
        dof_tables.emplace_back(*local_to_global_index_map_with_base_nodes_);
        dof_tables.emplace_back(*local_to_global_index_map_);
    }

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    GlobalExecutor::executeSelectedMemberDereferenced(
        global_assembler_, &VectorMatrixAssembler::assembleWithJacobian,
        local_assemblers_, pv.getActiveElementIDs(), dof_tables, t, dt, x, xdot,
        dxdot_dx, dx_dx, process_id, M, K, b, Jac, coupled_solutions_);

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
    if (use_monolithic_scheme_ || process_id == 0)
    {
        copyRhs(0, *hydraulic_flow_);
    }
    if (use_monolithic_scheme_ || process_id == 1)
    {
        copyRhs(1, *nodal_forces_);
    }
}

template <int DisplacementDim>
void HydroMechanicsProcess<DisplacementDim>::preTimestepConcreteProcess(
    std::vector<GlobalVector*> const& x, double const t, double const dt,
    const int process_id)
{
    DBUG("PreTimestep HydroMechanicsProcess.");

    if (hasMechanicalProcess(process_id))
    {
        ProcessLib::ProcessVariable const& pv =
            getProcessVariables(process_id)[0];
        GlobalExecutor::executeSelectedMemberOnDereferenced(
            &LocalAssemblerIF::preTimestep, local_assemblers_,
            pv.getActiveElementIDs(), *local_to_global_index_map_,
            *x[process_id], t, dt);
    }
}

template <int DisplacementDim>
void HydroMechanicsProcess<DisplacementDim>::postTimestepConcreteProcess(
    std::vector<GlobalVector*> const& x, double const t, double const dt,
    const int process_id)
{
    DBUG("PostTimestep HydroMechanicsProcess.");
    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];
    GlobalExecutor::executeSelectedMemberOnDereferenced(
        &LocalAssemblerIF::postTimestep, local_assemblers_,
        pv.getActiveElementIDs(), getDOFTable(process_id), *x[process_id], t,
        dt);
}

template <int DisplacementDim>
void HydroMechanicsProcess<DisplacementDim>::postNonLinearSolverConcreteProcess(
    GlobalVector const& x, const double t, double const dt,
    const int process_id)
{
    if (!hasMechanicalProcess(process_id))
    {
        return;
    }

    DBUG("PostNonLinearSolver HydroMechanicsProcess.");
    // Calculate strain, stress or other internal variables of mechanics.
    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];
    GlobalExecutor::executeSelectedMemberOnDereferenced(
        &LocalAssemblerIF::postNonLinearSolver, local_assemblers_,
        pv.getActiveElementIDs(), getDOFTable(process_id), x, t, dt,
        use_monolithic_scheme_);
}

template <int DisplacementDim>
void HydroMechanicsProcess<DisplacementDim>::computeSecondaryVariableConcrete(
    double const t, double const dt, GlobalVector const& x,
    GlobalVector const& x_dot, const int process_id)
{
    DBUG("Compute the secondary variables for HydroMechanicsProcess.");
    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];
    GlobalExecutor::executeSelectedMemberOnDereferenced(
        &LocalAssemblerIF::computeSecondaryVariable, local_assemblers_,
        pv.getActiveElementIDs(), getDOFTable(process_id), t, dt, x, x_dot,
        coupled_solutions_);
}

template <int DisplacementDim>
std::tuple<NumLib::LocalToGlobalIndexMap*, bool>
HydroMechanicsProcess<DisplacementDim>::getDOFTableForExtrapolatorData() const
{
    const bool manage_storage = false;
    return std::make_tuple(local_to_global_index_map_single_component_.get(),
                           manage_storage);
}

template <int DisplacementDim>
NumLib::LocalToGlobalIndexMap const&
HydroMechanicsProcess<DisplacementDim>::getDOFTable(const int process_id) const
{
    if (hasMechanicalProcess(process_id))
    {
        return *local_to_global_index_map_;
    }

    // For the equation of pressure
    return *local_to_global_index_map_with_base_nodes_;
}

template class HydroMechanicsProcess<2>;
template class HydroMechanicsProcess<3>;

}  // namespace HydroMechanics
}  // namespace ProcessLib
