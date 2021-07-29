/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
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
      _process_data(std::move(process_data))
{
    _nodal_forces = MeshLib::getOrCreateMeshProperty<double>(
        mesh, "NodalForces", MeshLib::MeshItemType::Node, DisplacementDim);

    _hydraulic_flow = MeshLib::getOrCreateMeshProperty<double>(
        mesh, "HydraulicFlow", MeshLib::MeshItemType::Node, 1);

    _integration_point_writer.emplace_back(
        std::make_unique<IntegrationPointWriter>(
            "sigma_ip",
            static_cast<int>(mesh.getDimension() == 2 ? 4 : 6) /*n components*/,
            integration_order,
            [this]()
            {
                // Result containing integration point data for each local
                // assembler.
                std::vector<std::vector<double>> result;
                result.resize(_local_assemblers.size());

                for (std::size_t i = 0; i < _local_assemblers.size(); ++i)
                {
                    auto const& local_asm = *_local_assemblers[i];

                    result[i] = local_asm.getSigma();
                }

                return result;
            }));
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
    const int mechanical_process_id = _use_monolithic_scheme ? 0 : 2;
    const int deformation_variable_id = _use_monolithic_scheme ? 2 : 0;
    ProcessLib::ThermoHydroMechanics::createLocalAssemblers<
        DisplacementDim, ThermoHydroMechanicsLocalAssembler>(
        mesh.getDimension(), mesh.getElements(), dof_table,
        // use displacement process variable to set shape function order
        getProcessVariables(mechanical_process_id)[deformation_variable_id]
            .get()
            .getShapeFunctionOrder(),
        _local_assemblers, mesh.isAxiallySymmetric(), integration_order,
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
        "velocity",
        makeExtrapolator(mesh.getDimension(), getExtrapolator(),
                         _local_assemblers,
                         &LocalAssemblerInterface::getIntPtDarcyVelocity));

    _process_data.pressure_interpolated =
        MeshLib::getOrCreateMeshProperty<double>(
            const_cast<MeshLib::Mesh&>(mesh), "pressure_interpolated",
            MeshLib::MeshItemType::Node, 1);

    _process_data.temperature_interpolated =
        MeshLib::getOrCreateMeshProperty<double>(
            const_cast<MeshLib::Mesh&>(mesh), "temperature_interpolated",
            MeshLib::MeshItemType::Node, 1);

    // Initialize local assemblers after all variables have been set.
    GlobalExecutor::executeMemberOnDereferenced(
        &LocalAssemblerInterface::initialize, _local_assemblers,
        *_local_to_global_index_map);

    // Set initial conditions for integration point data.
    for (auto const& ip_writer : _integration_point_writer)
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
        if (ip_meta_data.n_components !=
            mesh_property.getNumberOfGlobalComponents())
        {
            OGS_FATAL(
                "Different number of components in meta data ({:d}) than in "
                "the integration point field data for '{:s}': {:d}.",
                ip_meta_data.n_components, name,
                mesh_property.getNumberOfGlobalComponents());
        }

        // Now we have a properly named vtk's field data array and the
        // corresponding meta data.
        std::size_t position = 0;
        for (auto& local_asm : _local_assemblers)
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
}

template <int DisplacementDim>
void ThermoHydroMechanicsProcess<
    DisplacementDim>::initializeBoundaryConditions()
{
    if (_use_monolithic_scheme)
    {
        const int process_id_of_thermohydromechancs = 0;
        initializeProcessBoundaryConditionsAndSourceTerms(
            *_local_to_global_index_map, process_id_of_thermohydromechancs);
        return;
    }

    // Staggered scheme:
    // for the equations of heat transport
    const int thermal_process_id = 0;
    initializeProcessBoundaryConditionsAndSourceTerms(
        *_local_to_global_index_map_with_base_nodes, thermal_process_id);

    // for the equations of mass balance
    const int hydraulic_process_id = 1;
    initializeProcessBoundaryConditionsAndSourceTerms(
        *_local_to_global_index_map_with_base_nodes, hydraulic_process_id);

    // for the equations of deformation.
    const int mechanical_process_id = 2;
    initializeProcessBoundaryConditionsAndSourceTerms(
        *_local_to_global_index_map, mechanical_process_id);
}

template <int DisplacementDim>
void ThermoHydroMechanicsProcess<DisplacementDim>::assembleConcreteProcess(
    const double t, double const dt, std::vector<GlobalVector*> const& x,
    std::vector<GlobalVector*> const& xdot, int const process_id,
    GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b)
{
    DBUG("Assemble the equations for ThermoHydroMechanics");

    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_table = {std::ref(*_local_to_global_index_map)};

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    GlobalExecutor::executeSelectedMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assemble, _local_assemblers,
        pv.getActiveElementIDs(), dof_table, t, dt, x, xdot, process_id, M, K,
        b);
}

template <int DisplacementDim>
void ThermoHydroMechanicsProcess<DisplacementDim>::
    assembleWithJacobianConcreteProcess(
        const double t, double const dt, std::vector<GlobalVector*> const& x,
        std::vector<GlobalVector*> const& xdot, const double dxdot_dx,
        const double dx_dx, int const process_id, GlobalMatrix& M,
        GlobalMatrix& K, GlobalVector& b, GlobalMatrix& Jac)
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
        _local_assemblers, pv.getActiveElementIDs(), dof_tables, t, dt, x, xdot,
        dxdot_dx, dx_dx, process_id, M, K, b, Jac);

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
    if (_use_monolithic_scheme || process_id == 1)
    {
        copyRhs(0, *_hydraulic_flow);
    }
    if (_use_monolithic_scheme || process_id == 2)
    {
        copyRhs(1, *_nodal_forces);
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
        ProcessLib::ProcessVariable const& pv =
            getProcessVariables(process_id)[0];

        GlobalExecutor::executeSelectedMemberOnDereferenced(
            &LocalAssemblerInterface::preTimestep, _local_assemblers,
            pv.getActiveElementIDs(), *_local_to_global_index_map,
            *x[process_id], t, dt);
    }
}

template <int DisplacementDim>
void ThermoHydroMechanicsProcess<DisplacementDim>::postTimestepConcreteProcess(
    std::vector<GlobalVector*> const& x, double const t, double const dt,
    const int process_id)
{
    if (process_id != 0)
    {
        return;
    }

    DBUG("PostTimestep ThermoHydroMechanicsProcess.");
    std::vector<NumLib::LocalToGlobalIndexMap const*> dof_tables;
    auto const n_processes = x.size();
    dof_tables.reserve(n_processes);
    for (std::size_t process_id = 0; process_id < n_processes; ++process_id)
    {
        dof_tables.push_back(&getDOFTable(process_id));
    }

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    GlobalExecutor::executeSelectedMemberOnDereferenced(
        &LocalAssemblerInterface::postTimestep, _local_assemblers,
        pv.getActiveElementIDs(), dof_tables, x, t, dt);
}

template <int DisplacementDim>
void ThermoHydroMechanicsProcess<DisplacementDim>::
    postNonLinearSolverConcreteProcess(GlobalVector const& x,
                                       GlobalVector const& xdot, const double t,
                                       double const dt, const int process_id)
{
    if (!hasMechanicalProcess(process_id))
    {
        return;
    }

    DBUG("PostNonLinearSolver ThermoHydroMechanicsProcess.");
    // Calculate strain, stress or other internal variables of mechanics.

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    GlobalExecutor::executeSelectedMemberOnDereferenced(
        &LocalAssemblerInterface::postNonLinearSolver, _local_assemblers,
        pv.getActiveElementIDs(), getDOFTable(process_id), x, xdot, t, dt,
        _use_monolithic_scheme, process_id);
}

template <int DisplacementDim>
void ThermoHydroMechanicsProcess<DisplacementDim>::
    computeSecondaryVariableConcrete(double const t, double const dt,
                                     std::vector<GlobalVector*> const& x,
                                     GlobalVector const& x_dot,
                                     const int process_id)
{
    if (process_id != 0)
    {
        return;
    }

    DBUG("Compute the secondary variables for ThermoHydroMechanicsProcess.");
    std::vector<NumLib::LocalToGlobalIndexMap const*> dof_tables;
    auto const n_processes = x.size();
    dof_tables.reserve(n_processes);
    for (std::size_t process_id = 0; process_id < n_processes; ++process_id)
    {
        dof_tables.push_back(&getDOFTable(process_id));
    }

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];
    GlobalExecutor::executeSelectedMemberOnDereferenced(
        &LocalAssemblerInterface::computeSecondaryVariable, _local_assemblers,
        pv.getActiveElementIDs(), dof_tables, t, dt, x, x_dot, process_id);
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
