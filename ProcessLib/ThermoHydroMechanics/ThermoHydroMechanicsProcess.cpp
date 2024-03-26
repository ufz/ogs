/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ThermoHydroMechanicsProcess.h"

#include <cassert>

#include "MeshLib/Elements/Utils.h"
#include "MeshLib/Utils/getOrCreateMeshProperty.h"
#include "NumLib/DOF/ComputeSparsityPattern.h"
#include "NumLib/DOF/DOFTableUtil.h"
#include "ProcessLib/Deformation/SolidMaterialInternalToSecondaryVariables.h"
#include "ProcessLib/Process.h"
#include "ProcessLib/Utils/CreateLocalAssemblersTaylorHood.h"
#include "ProcessLib/Utils/SetIPDataInitialConditions.h"
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
      _process_data(std::move(process_data))
{
    _nodal_forces = MeshLib::getOrCreateMeshProperty<double>(
        mesh, "NodalForces", MeshLib::MeshItemType::Node, DisplacementDim);

    _hydraulic_flow = MeshLib::getOrCreateMeshProperty<double>(
        mesh, "MassFlowRate", MeshLib::MeshItemType::Node, 1);

    _heat_flux = MeshLib::getOrCreateMeshProperty<double>(
        mesh, "HeatFlowRate", MeshLib::MeshItemType::Node, 1);

    _integration_point_writer.emplace_back(
        std::make_unique<MeshLib::IntegrationPointWriter>(
            "sigma_ip",
            static_cast<int>(mesh.getDimension() == 2 ? 4 : 6) /*n components*/,
            integration_order, _local_assemblers,
            &LocalAssemblerInterface<DisplacementDim>::getSigma));

    _integration_point_writer.emplace_back(
        std::make_unique<MeshLib::IntegrationPointWriter>(
            "epsilon_m_ip",
            static_cast<int>(mesh.getDimension() == 2 ? 4 : 6) /*n components*/,
            integration_order, _local_assemblers,
            &LocalAssemblerInterface<DisplacementDim>::getEpsilonM));

    _integration_point_writer.emplace_back(
        std::make_unique<MeshLib::IntegrationPointWriter>(
            "epsilon_ip",
            static_cast<int>(mesh.getDimension() == 2 ? 4 : 6) /*n components*/,
            integration_order, _local_assemblers,
            &LocalAssemblerInterface<DisplacementDim>::getEpsilon));
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
    if (_use_monolithic_scheme || process_id == 2)
    {
        auto const& l = *_local_to_global_index_map;
        return {l.dofSizeWithoutGhosts(), l.dofSizeWithoutGhosts(),
                &l.getGhostIndices(), &this->_sparsity_pattern};
    }

    // For staggered scheme and T or H process (pressure).
    auto const& l = *_local_to_global_index_map_with_base_nodes;
    return {l.dofSizeWithoutGhosts(), l.dofSizeWithoutGhosts(),
            &l.getGhostIndices(), &_sparsity_pattern_with_linear_element};
}

template <int DisplacementDim>
void ThermoHydroMechanicsProcess<DisplacementDim>::constructDofTable()
{
    // Create single component dof in every of the mesh's nodes.
    _mesh_subset_all_nodes =
        std::make_unique<MeshLib::MeshSubset>(_mesh, _mesh.getNodes());
    // Create single component dof in the mesh's base nodes.
    _base_nodes = MeshLib::getBaseNodes(_mesh.getElements());
    _mesh_subset_base_nodes =
        std::make_unique<MeshLib::MeshSubset>(_mesh, _base_nodes);

    // TODO move the two data members somewhere else.
    // for extrapolation of secondary variables of stress or strain
    std::vector<MeshLib::MeshSubset> all_mesh_subsets_single_component{
        *_mesh_subset_all_nodes};
    _local_to_global_index_map_single_component =
        std::make_unique<NumLib::LocalToGlobalIndexMap>(
            std::move(all_mesh_subsets_single_component),
            // by location order is needed for output
            NumLib::ComponentOrder::BY_LOCATION);

    if (_use_monolithic_scheme)
    {
        // For temperature, which is the first
        std::vector<MeshLib::MeshSubset> all_mesh_subsets{
            *_mesh_subset_base_nodes};

        // For pressure, which is the second
        all_mesh_subsets.push_back(*_mesh_subset_base_nodes);

        // For displacement.
        const int monolithic_process_id = 0;
        std::generate_n(std::back_inserter(all_mesh_subsets),
                        getProcessVariables(monolithic_process_id)[2]
                            .get()
                            .getNumberOfGlobalComponents(),
                        [&]() { return *_mesh_subset_all_nodes; });

        std::vector<int> const vec_n_components{1, 1, DisplacementDim};
        _local_to_global_index_map =
            std::make_unique<NumLib::LocalToGlobalIndexMap>(
                std::move(all_mesh_subsets), vec_n_components,
                NumLib::ComponentOrder::BY_LOCATION);
        assert(_local_to_global_index_map);
    }
    else
    {
        // For displacement equation.
        const int process_id = 2;
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

        // For pressure equation or temperature equation.
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
void ThermoHydroMechanicsProcess<DisplacementDim>::initializeConcreteProcess(
    NumLib::LocalToGlobalIndexMap const& dof_table,
    MeshLib::Mesh const& mesh,
    unsigned const integration_order)
{
    ProcessLib::createLocalAssemblersHM<DisplacementDim,
                                        ThermoHydroMechanicsLocalAssembler>(
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

    add_secondary_variable(
        "sigma",
        MathLib::KelvinVector::KelvinVectorType<
            DisplacementDim>::RowsAtCompileTime,
        &LocalAssemblerInterface<DisplacementDim>::getIntPtSigma);

    add_secondary_variable(
        "sigma_ice",
        MathLib::KelvinVector::KelvinVectorType<
            DisplacementDim>::RowsAtCompileTime,
        &LocalAssemblerInterface<DisplacementDim>::getIntPtSigmaIce);

    add_secondary_variable(
        "epsilon_m",
        MathLib::KelvinVector::KelvinVectorType<
            DisplacementDim>::RowsAtCompileTime,
        &LocalAssemblerInterface<DisplacementDim>::getIntPtEpsilonM);

    add_secondary_variable(
        "epsilon",
        MathLib::KelvinVector::KelvinVectorType<
            DisplacementDim>::RowsAtCompileTime,
        &LocalAssemblerInterface<DisplacementDim>::getIntPtEpsilon);

    add_secondary_variable(
        "ice_volume_fraction", 1,
        &LocalAssemblerInterface<DisplacementDim>::getIntPtIceVolume);

    add_secondary_variable(
        "velocity", mesh.getDimension(),
        &LocalAssemblerInterface<DisplacementDim>::getIntPtDarcyVelocity);

    add_secondary_variable(
        "fluid_density", 1,
        &LocalAssemblerInterface<DisplacementDim>::getIntPtFluidDensity);

    add_secondary_variable(
        "viscosity", 1,
        &LocalAssemblerInterface<DisplacementDim>::getIntPtViscosity);

    _process_data.element_fluid_density =
        MeshLib::getOrCreateMeshProperty<double>(
            const_cast<MeshLib::Mesh&>(mesh), "fluid_density_avg",
            MeshLib::MeshItemType::Cell, 1);

    _process_data.element_viscosity = MeshLib::getOrCreateMeshProperty<double>(
        const_cast<MeshLib::Mesh&>(mesh), "viscosity_avg",
        MeshLib::MeshItemType::Cell, 1);

    _process_data.element_stresses = MeshLib::getOrCreateMeshProperty<double>(
        const_cast<MeshLib::Mesh&>(mesh), "stress_avg",
        MeshLib::MeshItemType::Cell,
        MathLib::KelvinVector::KelvinVectorType<
            DisplacementDim>::RowsAtCompileTime);

    _process_data.pressure_interpolated =
        MeshLib::getOrCreateMeshProperty<double>(
            const_cast<MeshLib::Mesh&>(mesh), "pressure_interpolated",
            MeshLib::MeshItemType::Node, 1);

    _process_data.temperature_interpolated =
        MeshLib::getOrCreateMeshProperty<double>(
            const_cast<MeshLib::Mesh&>(mesh), "temperature_interpolated",
            MeshLib::MeshItemType::Node, 1);

    //
    // enable output of internal variables defined by material models
    //
    ProcessLib::Deformation::solidMaterialInternalToSecondaryVariables<
        LocalAssemblerInterface<DisplacementDim>>(_process_data.solid_materials,
                                                  add_secondary_variable);

    ProcessLib::Deformation::
        solidMaterialInternalVariablesToIntegrationPointWriter(
            _process_data.solid_materials, _local_assemblers,
            _integration_point_writer, integration_order);

    setIPDataInitialConditions(_integration_point_writer, mesh.getProperties(),
                               _local_assemblers);

    // Initialize local assemblers after all variables have been set.
    GlobalExecutor::executeMemberOnDereferenced(
        &LocalAssemblerInterface<DisplacementDim>::initialize,
        _local_assemblers, *_local_to_global_index_map);
}

template <int DisplacementDim>
void ThermoHydroMechanicsProcess<DisplacementDim>::initializeBoundaryConditions(
    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const& media)
{
    if (_use_monolithic_scheme)
    {
        const int process_id_of_thermohydromechancs = 0;
        initializeProcessBoundaryConditionsAndSourceTerms(
            *_local_to_global_index_map, process_id_of_thermohydromechancs,
            media);
        return;
    }

    // Staggered scheme:
    // for the equations of heat transport
    const int thermal_process_id = 0;
    initializeProcessBoundaryConditionsAndSourceTerms(
        *_local_to_global_index_map_with_base_nodes, thermal_process_id, media);

    // for the equations of mass balance
    const int hydraulic_process_id = 1;
    initializeProcessBoundaryConditionsAndSourceTerms(
        *_local_to_global_index_map_with_base_nodes, hydraulic_process_id,
        media);

    // for the equations of deformation.
    const int mechanical_process_id = 2;
    initializeProcessBoundaryConditionsAndSourceTerms(
        *_local_to_global_index_map, mechanical_process_id, media);
}

template <int DisplacementDim>
void ThermoHydroMechanicsProcess<DisplacementDim>::
    setInitialConditionsConcreteProcess(std::vector<GlobalVector*>& x,
                                        double const t,
                                        int const process_id)
{
    DBUG("SetInitialConditions ThermoHydroMechanicsProcess.");

    auto get_a_dof_table_func = [this](const int process_id) -> auto&
    {
        return getDOFTable(process_id);
    };
    GlobalExecutor::executeMemberOnDereferenced(
        &LocalAssemblerInterface<DisplacementDim>::setInitialConditions,
        _local_assemblers, NumLib::getDOFTables(x.size(), get_a_dof_table_func),
        x, t, process_id);
}

template <int DisplacementDim>
void ThermoHydroMechanicsProcess<DisplacementDim>::assembleConcreteProcess(
    const double t, double const dt, std::vector<GlobalVector*> const& x,
    std::vector<GlobalVector*> const& x_prev, int const process_id,
    GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b)
{
    DBUG("Assemble the equations for ThermoHydroMechanics");

    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_table = {std::ref(*_local_to_global_index_map)};

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    GlobalExecutor::executeSelectedMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assemble, _local_assemblers,
        pv.getActiveElementIDs(), dof_table, t, dt, x, x_prev, process_id, M, K,
        b);
}

template <int DisplacementDim>
void ThermoHydroMechanicsProcess<DisplacementDim>::
    assembleWithJacobianConcreteProcess(
        const double t, double const dt, std::vector<GlobalVector*> const& x,
        std::vector<GlobalVector*> const& x_prev, int const process_id,
        GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b, GlobalMatrix& Jac)
{
    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_tables;
    // For the monolithic scheme
    if (_use_monolithic_scheme)
    {
        DBUG(
            "Assemble the Jacobian of ThermoHydroMechanics for the monolithic "
            "scheme.");
        dof_tables.emplace_back(*_local_to_global_index_map);
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
        dof_tables.emplace_back(*_local_to_global_index_map_with_base_nodes);
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
        if (_use_monolithic_scheme)
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
    if (_use_monolithic_scheme || process_id == 0)
    {
        copyRhs(0, *_heat_flux);
    }
    if (_use_monolithic_scheme || process_id == 1)
    {
        copyRhs(1, *_hydraulic_flow);
    }
    if (_use_monolithic_scheme || process_id == 2)
    {
        copyRhs(2, *_nodal_forces);
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
            &LocalAssemblerInterface<DisplacementDim>::preTimestep,
            _local_assemblers, *_local_to_global_index_map, *x[process_id], t,
            dt);
    }
}

template <int DisplacementDim>
void ThermoHydroMechanicsProcess<DisplacementDim>::postTimestepConcreteProcess(
    std::vector<GlobalVector*> const& x,
    std::vector<GlobalVector*> const& x_prev, double const t, double const dt,
    const int process_id)
{
    if (process_id != 0)
    {
        return;
    }

    DBUG("PostTimestep ThermoHydroMechanicsProcess.");

    auto get_a_dof_table_func = [this](const int processe_id) -> auto&
    {
        return getDOFTable(processe_id);
    };
    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];
    GlobalExecutor::executeSelectedMemberOnDereferenced(
        &LocalAssemblerInterface<DisplacementDim>::postTimestep,
        _local_assemblers, pv.getActiveElementIDs(),
        NumLib::getDOFTables(x.size(), get_a_dof_table_func), x, x_prev, t, dt,
        process_id);
}

template <int DisplacementDim>
void ThermoHydroMechanicsProcess<DisplacementDim>::
    computeSecondaryVariableConcrete(double const t, double const dt,
                                     std::vector<GlobalVector*> const& x,
                                     GlobalVector const& x_prev,
                                     const int process_id)
{
    if (process_id != 0)
    {
        return;
    }

    DBUG("Compute the secondary variables for ThermoHydroMechanicsProcess.");

    auto get_a_dof_table_func = [this](const int processe_id) -> auto&
    {
        return getDOFTable(processe_id);
    };
    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];
    GlobalExecutor::executeSelectedMemberOnDereferenced(
        &LocalAssemblerInterface<DisplacementDim>::computeSecondaryVariable,
        _local_assemblers, pv.getActiveElementIDs(),
        NumLib::getDOFTables(x.size(), get_a_dof_table_func), t, dt, x, x_prev,
        process_id);
}

template <int DisplacementDim>
std::tuple<NumLib::LocalToGlobalIndexMap*, bool> ThermoHydroMechanicsProcess<
    DisplacementDim>::getDOFTableForExtrapolatorData() const
{
    const bool manage_storage = false;
    return std::make_tuple(_local_to_global_index_map_single_component.get(),
                           manage_storage);
}

template <int DisplacementDim>
NumLib::LocalToGlobalIndexMap const&
ThermoHydroMechanicsProcess<DisplacementDim>::getDOFTable(
    const int process_id) const
{
    if (hasMechanicalProcess(process_id))
    {
        return *_local_to_global_index_map;
    }

    // For the equation of pressure
    return *_local_to_global_index_map_with_base_nodes;
}

template class ThermoHydroMechanicsProcess<2>;
template class ThermoHydroMechanicsProcess<3>;

}  // namespace ThermoHydroMechanics
}  // namespace ProcessLib
