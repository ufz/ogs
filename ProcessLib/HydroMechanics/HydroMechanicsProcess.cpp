/**
 * \file
 * \copyright
 * Copyright (c) 2012-2023, OpenGeoSys Community (http://www.opengeosys.org)
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
#include "MeshLib/Utils/getOrCreateMeshProperty.h"
#include "NumLib/DOF/ComputeSparsityPattern.h"
#include "ProcessLib/Deformation/SolidMaterialInternalToSecondaryVariables.h"
#include "ProcessLib/Process.h"
#include "ProcessLib/Utils/CreateLocalAssemblersTaylorHood.h"
#include "ProcessLib/Utils/SetIPDataInitialConditions.h"

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
      _process_data(std::move(process_data))
{
    _nodal_forces = MeshLib::getOrCreateMeshProperty<double>(
        mesh, "NodalForces", MeshLib::MeshItemType::Node, DisplacementDim);

    _hydraulic_flow = MeshLib::getOrCreateMeshProperty<double>(
        mesh, "MassFlowRate", MeshLib::MeshItemType::Node, 1);

    _integration_point_writer.emplace_back(
        std::make_unique<MeshLib::IntegrationPointWriter>(
            "sigma_ip",
            static_cast<int>(mesh.getDimension() == 2 ? 4 : 6) /*n components*/,
            integration_order, _local_assemblers, &LocalAssemblerIF::getSigma));

    _integration_point_writer.emplace_back(
        std::make_unique<MeshLib::IntegrationPointWriter>(
            "epsilon_ip",
            static_cast<int>(mesh.getDimension() == 2 ? 4 : 6) /*n components*/,
            integration_order, _local_assemblers,
            &LocalAssemblerIF::getEpsilon));

    if (!_use_monolithic_scheme)
    {
        _integration_point_writer.emplace_back(
            std::make_unique<MeshLib::IntegrationPointWriter>(
                "strain_rate_variable_ip", 1, integration_order,
                _local_assemblers, &LocalAssemblerIF::getStrainRateVariable));
    }
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
    if (process_id == _process_data.mechanics_related_process_id)
    {
        auto const& l = *_local_to_global_index_map;
        return {l.dofSizeWithoutGhosts(), l.dofSizeWithoutGhosts(),
                &l.getGhostIndices(), &this->_sparsity_pattern};
    }

    // For staggered scheme and H process (pressure).
    auto const& l = *_local_to_global_index_map_with_base_nodes;
    return {l.dofSizeWithoutGhosts(), l.dofSizeWithoutGhosts(),
            &l.getGhostIndices(), &_sparsity_pattern_with_linear_element};
}

template <int DisplacementDim>
void HydroMechanicsProcess<DisplacementDim>::constructDofTable()
{
    // Create single component dof in every of the mesh's nodes.
    _mesh_subset_all_nodes = std::make_unique<MeshLib::MeshSubset>(
        _mesh, _mesh.getNodes(), _process_data.use_taylor_hood_elements);

    // Create single component dof in the mesh's base nodes.
    _base_nodes = MeshLib::getBaseNodes(_mesh.getElements());
    _mesh_subset_base_nodes = std::make_unique<MeshLib::MeshSubset>(
        _mesh, _base_nodes, _process_data.use_taylor_hood_elements);

    // TODO move the two data members somewhere else.
    // for extrapolation of secondary variables of stress or strain
    std::vector<MeshLib::MeshSubset> all_mesh_subsets_single_component{
        *_mesh_subset_all_nodes};
    _local_to_global_index_map_single_component =
        std::make_unique<NumLib::LocalToGlobalIndexMap>(
            std::move(all_mesh_subsets_single_component),
            // by location order is needed for output
            NumLib::ComponentOrder::BY_LOCATION);

    if (_process_data.isMonolithicSchemeUsed())
    {
        // For pressure, which is the first
        std::vector<MeshLib::MeshSubset> all_mesh_subsets{
            *_mesh_subset_base_nodes};

        // For displacement.
        const int monolithic_process_id = 0;
        std::generate_n(std::back_inserter(all_mesh_subsets),
                        getProcessVariables(monolithic_process_id)[1]
                            .get()
                            .getNumberOfGlobalComponents(),
                        [&]() { return *_mesh_subset_all_nodes; });

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
        std::vector<MeshLib::MeshSubset> all_mesh_subsets;
        std::generate_n(std::back_inserter(all_mesh_subsets),
                        getProcessVariables(process_id)[0]
                            .get()
                            .getNumberOfGlobalComponents(),
                        [&]() { return *_mesh_subset_all_nodes; });

        std::vector<int> const vec_n_components{DisplacementDim};
        _local_to_global_index_map =
            std::make_unique<NumLib::LocalToGlobalIndexMap>(
                std::move(all_mesh_subsets), vec_n_components,
                NumLib::ComponentOrder::BY_LOCATION);

        // For pressure equation.
        // Collect the mesh subsets with base nodes in a vector.
        std::vector<MeshLib::MeshSubset> all_mesh_subsets_base_nodes{
            *_mesh_subset_base_nodes};
        _local_to_global_index_map_with_base_nodes =
            std::make_unique<NumLib::LocalToGlobalIndexMap>(
                std::move(all_mesh_subsets_base_nodes),
                // by location order is needed for output
                NumLib::ComponentOrder::BY_LOCATION);

        _sparsity_pattern_with_linear_element = NumLib::computeSparsityPattern(
            *_local_to_global_index_map_with_base_nodes, _mesh);

        assert(_local_to_global_index_map);
        assert(_local_to_global_index_map_with_base_nodes);
    }
}

template <int DisplacementDim>
void HydroMechanicsProcess<DisplacementDim>::initializeConcreteProcess(
    NumLib::LocalToGlobalIndexMap const& dof_table,
    MeshLib::Mesh const& mesh,
    unsigned const integration_order)
{
    ProcessLib::createLocalAssemblersHM<DisplacementDim,
                                        HydroMechanicsLocalAssembler>(
        mesh.getElements(), dof_table, _local_assemblers,
        NumLib::IntegrationOrder{integration_order}, mesh.isAxiallySymmetric(),
        _process_data);

    auto add_secondary_variable = [&](std::string const& name,
                                      int const num_components,
                                      auto get_ip_values_function)
    {
        _secondary_variables.addSecondaryVariable(
            name,
            makeExtrapolator(num_components, getExtrapolator(),
                             _local_assemblers,
                             std::move(get_ip_values_function)));
    };

    add_secondary_variable("sigma",
                           MathLib::KelvinVector::KelvinVectorType<
                               DisplacementDim>::RowsAtCompileTime,
                           &LocalAssemblerIF::getIntPtSigma);

    add_secondary_variable("epsilon",
                           MathLib::KelvinVector::KelvinVectorType<
                               DisplacementDim>::RowsAtCompileTime,
                           &LocalAssemblerIF::getIntPtEpsilon);

    add_secondary_variable("velocity", DisplacementDim,
                           &LocalAssemblerIF::getIntPtDarcyVelocity);

    //
    // enable output of internal variables defined by material models
    //
    ProcessLib::Deformation::solidMaterialInternalToSecondaryVariables<
        LocalAssemblerIF>(_process_data.solid_materials,
                          add_secondary_variable);

    _process_data.pressure_interpolated =
        MeshLib::getOrCreateMeshProperty<double>(
            const_cast<MeshLib::Mesh&>(mesh), "pressure_interpolated",
            MeshLib::MeshItemType::Node, 1);

    _process_data.principal_stress_vector[0] =
        MeshLib::getOrCreateMeshProperty<double>(
            const_cast<MeshLib::Mesh&>(mesh), "principal_stress_vector_1",
            MeshLib::MeshItemType::Cell, 3);

    _process_data.principal_stress_vector[1] =
        MeshLib::getOrCreateMeshProperty<double>(
            const_cast<MeshLib::Mesh&>(mesh), "principal_stress_vector_2",
            MeshLib::MeshItemType::Cell, 3);

    _process_data.principal_stress_vector[2] =
        MeshLib::getOrCreateMeshProperty<double>(
            const_cast<MeshLib::Mesh&>(mesh), "principal_stress_vector_3",
            MeshLib::MeshItemType::Cell, 3);

    _process_data.principal_stress_values =
        MeshLib::getOrCreateMeshProperty<double>(
            const_cast<MeshLib::Mesh&>(mesh), "principal_stress_values",
            MeshLib::MeshItemType::Cell, 3);

    _process_data.permeability = MeshLib::getOrCreateMeshProperty<double>(
        const_cast<MeshLib::Mesh&>(mesh), "permeability",
        MeshLib::MeshItemType::Cell,
        MathLib::KelvinVector::kelvin_vector_dimensions(DisplacementDim));

    setIPDataInitialConditions(_integration_point_writer, mesh.getProperties(),
                               _local_assemblers);

    // Initialize local assemblers after all variables have been set.
    GlobalExecutor::executeMemberOnDereferenced(&LocalAssemblerIF::initialize,
                                                _local_assemblers,
                                                *_local_to_global_index_map);
}

template <int DisplacementDim>
void HydroMechanicsProcess<DisplacementDim>::initializeBoundaryConditions(
    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const& media)
{
    if (_process_data.isMonolithicSchemeUsed())
    {
        const int process_id_of_hydromechanics = 0;
        initializeProcessBoundaryConditionsAndSourceTerms(
            *_local_to_global_index_map, process_id_of_hydromechanics, media);
        return;
    }

    // Staggered scheme:
    // for the equations of pressure
    const int hydraulic_process_id = 0;
    initializeProcessBoundaryConditionsAndSourceTerms(
        *_local_to_global_index_map_with_base_nodes, hydraulic_process_id,
        media);

    // for the equations of deformation.
    const int mechanical_process_id = 1;
    initializeProcessBoundaryConditionsAndSourceTerms(
        *_local_to_global_index_map, mechanical_process_id, media);
}

template <int DisplacementDim>
void HydroMechanicsProcess<DisplacementDim>::assembleConcreteProcess(
    const double t, double const dt, std::vector<GlobalVector*> const& x,
    std::vector<GlobalVector*> const& x_prev, int const process_id,
    GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b)
{
    DBUG("Assemble the equations for HydroMechanics");

    // Note: This assembly function is for the Picard nonlinear solver. Since
    // only the Newton-Raphson method is employed to simulate coupled HM
    // processes in this class, this function is actually not used so far.

    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_table = {std::ref(*_local_to_global_index_map)};

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    GlobalExecutor::executeSelectedMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assemble, _local_assemblers,
        pv.getActiveElementIDs(), dof_table, t, dt, x, x_prev, process_id, M, K,
        b);
}

template <int DisplacementDim>
void HydroMechanicsProcess<DisplacementDim>::
    assembleWithJacobianConcreteProcess(
        const double t, double const dt, std::vector<GlobalVector*> const& x,
        std::vector<GlobalVector*> const& x_prev, int const process_id,
        GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b, GlobalMatrix& Jac)
{
    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_tables;
    // For the monolithic scheme
    bool const use_monolithic_scheme = _process_data.isMonolithicSchemeUsed();
    if (use_monolithic_scheme)
    {
        DBUG(
            "Assemble the Jacobian of HydroMechanics for the monolithic "
            "scheme.");
        dof_tables.emplace_back(*_local_to_global_index_map);
    }
    else
    {
        // For the staggered scheme
        if (process_id == _process_data.hydraulic_process_id)
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
        dof_tables.emplace_back(*_local_to_global_index_map_with_base_nodes);
        dof_tables.emplace_back(*_local_to_global_index_map);
    }

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    GlobalExecutor::executeSelectedMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assembleWithJacobian,
        _local_assemblers, pv.getActiveElementIDs(), dof_tables, t, dt, x,
        x_prev, process_id, M, K, b, Jac);

    auto copyRhs = [&](int const variable_id, auto& output_vector)
    {
        if (use_monolithic_scheme)
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
    if (process_id == _process_data.hydraulic_process_id)
    {
        copyRhs(0, *_hydraulic_flow);
    }
    if (process_id == _process_data.mechanics_related_process_id)
    {
        copyRhs(1, *_nodal_forces);
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
            &LocalAssemblerIF::preTimestep, _local_assemblers,
            pv.getActiveElementIDs(), *_local_to_global_index_map,
            *x[process_id], t, dt);
    }
}

template <int DisplacementDim>
void HydroMechanicsProcess<DisplacementDim>::postTimestepConcreteProcess(
    std::vector<GlobalVector*> const& x,
    std::vector<GlobalVector*> const& x_prev, double const t, double const dt,
    const int process_id)
{
    if (process_id != _process_data.hydraulic_process_id)
    {
        return;
    }

    DBUG("PostTimestep HydroMechanicsProcess.");
    std::vector<NumLib::LocalToGlobalIndexMap const*> dof_tables;
    auto const n_processes = x.size();
    dof_tables.reserve(n_processes);
    for (std::size_t process_id = 0; process_id < n_processes; ++process_id)
    {
        dof_tables.push_back(&getDOFTable(process_id));
    }

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];
    GlobalExecutor::executeSelectedMemberOnDereferenced(
        &LocalAssemblerIF::postTimestep, _local_assemblers,
        pv.getActiveElementIDs(), dof_tables, x, x_prev, t, dt,
        false /*unused*/, process_id);
}

template <int DisplacementDim>
void HydroMechanicsProcess<DisplacementDim>::postNonLinearSolverConcreteProcess(
    std::vector<GlobalVector*> const& x,
    std::vector<GlobalVector*> const& x_prev, const double t, double const dt,
    const int process_id)
{
    DBUG("PostNonLinearSolver HydroMechanicsProcess.");

    // Calculate strain, stress or other internal variables of mechanics.
    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];
    GlobalExecutor::executeSelectedMemberOnDereferenced(
        &LocalAssemblerIF::postNonLinearSolver, _local_assemblers,
        pv.getActiveElementIDs(), getDOFTable(process_id), x, x_prev, t, dt,
        false /*unused*/, process_id);
}

template <int DisplacementDim>
void HydroMechanicsProcess<DisplacementDim>::
    setInitialConditionsConcreteProcess(std::vector<GlobalVector*>& x,
                                        double const t,
                                        int const process_id)
{
    if (process_id != _process_data.hydraulic_process_id)
    {
        return;
    }

    DBUG("Set initial conditions of HydroMechanicsProcess.");

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];
    GlobalExecutor::executeSelectedMemberOnDereferenced(
        &LocalAssemblerIF::setInitialConditions, _local_assemblers,
        pv.getActiveElementIDs(), getDOFTable(process_id), *x[process_id], t,
        _process_data.isMonolithicSchemeUsed(), process_id);
}

template <int DisplacementDim>
void HydroMechanicsProcess<DisplacementDim>::computeSecondaryVariableConcrete(
    double const t, double const dt, std::vector<GlobalVector*> const& x,
    GlobalVector const& x_prev, const int process_id)
{
    if (process_id != _process_data.hydraulic_process_id)
    {
        return;
    }

    DBUG("Compute the secondary variables for HydroMechanicsProcess.");
    std::vector<NumLib::LocalToGlobalIndexMap const*> dof_tables;
    auto const n_processes = x.size();
    dof_tables.reserve(n_processes);
    for (std::size_t process_id = 0; process_id < n_processes; ++process_id)
    {
        dof_tables.push_back(&getDOFTable(process_id));
    }

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];
    GlobalExecutor::executeSelectedMemberOnDereferenced(
        &LocalAssemblerIF::computeSecondaryVariable, _local_assemblers,
        pv.getActiveElementIDs(), dof_tables, t, dt, x, x_prev, process_id);
}

template <int DisplacementDim>
std::tuple<NumLib::LocalToGlobalIndexMap*, bool>
HydroMechanicsProcess<DisplacementDim>::getDOFTableForExtrapolatorData() const
{
    const bool manage_storage = false;
    return std::make_tuple(_local_to_global_index_map_single_component.get(),
                           manage_storage);
}

template <int DisplacementDim>
NumLib::LocalToGlobalIndexMap const&
HydroMechanicsProcess<DisplacementDim>::getDOFTable(const int process_id) const
{
    if (hasMechanicalProcess(process_id))
    {
        return *_local_to_global_index_map;
    }

    // For the equation of pressure
    return *_local_to_global_index_map_with_base_nodes;
}

template class HydroMechanicsProcess<2>;
template class HydroMechanicsProcess<3>;

}  // namespace HydroMechanics
}  // namespace ProcessLib
