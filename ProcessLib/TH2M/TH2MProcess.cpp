/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "TH2MProcess.h"

#include <cassert>

#include "MeshLib/Elements/Utils.h"
#include "NumLib/DOF/ComputeSparsityPattern.h"
#include "ProcessLib/Process.h"
#include "ProcessLib/TH2M/CreateLocalAssemblers.h"
#include "TH2MFEM.h"
#include "TH2MProcessData.h"

namespace ProcessLib
{
namespace TH2M
{
template <int DisplacementDim>
TH2MProcess<DisplacementDim>::TH2MProcess(
    std::string name,
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    unsigned const integration_order,
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
        process_variables,
    TH2MProcessData<DisplacementDim>&& process_data,
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
}

template <int DisplacementDim>
bool TH2MProcess<DisplacementDim>::isLinear() const
{
    return false;
}

template <int DisplacementDim>
MathLib::MatrixSpecifications
TH2MProcess<DisplacementDim>::getMatrixSpecifications(
    const int process_id) const
{
    // For the monolithic scheme or the M process (deformation) in the staggered
    // scheme.
    if (_use_monolithic_scheme || process_id == deformation_process_id)
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
void TH2MProcess<DisplacementDim>::constructDofTable()
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
        // For gas pressure, which is the first
        std::vector<MeshLib::MeshSubset> all_mesh_subsets{
            *_mesh_subset_base_nodes};

        // For capillary pressure, which is the second
        all_mesh_subsets.push_back(*_mesh_subset_base_nodes);

        // For temperature, which is the third
        all_mesh_subsets.push_back(*_mesh_subset_base_nodes);

        // For displacement.
        std::generate_n(
            std::back_inserter(all_mesh_subsets),
            getProcessVariables(monolithic_process_id)[deformation_process_id]
                .get()
                .getNumberOfGlobalComponents(),
            [&]() { return *_mesh_subset_all_nodes; });

        std::vector<int> const vec_n_components{
            n_gas_pressure_components, n_capillary_pressure_components,
            n_temperature_components, n_displacement_components};

        _local_to_global_index_map =
            std::make_unique<NumLib::LocalToGlobalIndexMap>(
                std::move(all_mesh_subsets), vec_n_components,
                NumLib::ComponentOrder::BY_LOCATION);
        assert(_local_to_global_index_map);
    }
    else
    {
        OGS_FATAL("A Staggered version of TH2M is not implemented.");
    }
}

template <int DisplacementDim>
void TH2MProcess<DisplacementDim>::initializeConcreteProcess(
    NumLib::LocalToGlobalIndexMap const& dof_table,
    MeshLib::Mesh const& mesh,
    unsigned const integration_order)
{
    const int mechanical_process_id =
        _use_monolithic_scheme ? 0 : deformation_process_id;
    const int deformation_variable_id =
        _use_monolithic_scheme ? deformation_process_id : 0;
    ProcessLib::TH2M::createLocalAssemblers<DisplacementDim,
                                            TH2MLocalAssembler>(
        mesh.getDimension(), mesh.getElements(), dof_table,
        // use displacement process variable to set shape function order
        getProcessVariables(mechanical_process_id)[deformation_variable_id]
            .get()
            .getShapeFunctionOrder(),
        _local_assemblers, mesh.isAxiallySymmetric(), integration_order,
        _process_data);

    auto add_secondary_variable = [&](std::string const& name,
                                      int const num_components,
                                      auto get_ip_values_function) {
        _secondary_variables.addSecondaryVariable(
            name,
            makeExtrapolator(num_components, getExtrapolator(),
                             _local_assemblers,
                             std::move(get_ip_values_function)));
    };

    add_secondary_variable("sigma",
                           MathLib::KelvinVector::KelvinVectorType<
                               DisplacementDim>::RowsAtCompileTime,
                           &LocalAssemblerInterface::getIntPtSigma);
    add_secondary_variable("epsilon",
                           MathLib::KelvinVector::KelvinVectorType<
                               DisplacementDim>::RowsAtCompileTime,
                           &LocalAssemblerInterface::getIntPtEpsilon);
    add_secondary_variable("velocity_gas", mesh.getDimension(),
                           &LocalAssemblerInterface::getIntPtDarcyVelocityGas);
    add_secondary_variable(
        "velocity_liquid", mesh.getDimension(),
        &LocalAssemblerInterface::getIntPtDarcyVelocityLiquid);
    add_secondary_variable("saturation", 1,
                           &LocalAssemblerInterface::getIntPtSaturation);
    add_secondary_variable("vapour_pressure", 1,
                           &LocalAssemblerInterface::getIntPtVapourPressure);
    add_secondary_variable("porosity", 1,
                           &LocalAssemblerInterface::getIntPtPorosity);
    add_secondary_variable("gas_density", 1,
                           &LocalAssemblerInterface::getIntPtGasDensity);
    add_secondary_variable("solid_density", 1,
                           &LocalAssemblerInterface::getIntPtSolidDensity);
    add_secondary_variable("liquid_density", 1,
                           &LocalAssemblerInterface::getIntPtLiquidDensity);
    add_secondary_variable("liquid_pressure", 1,
                           &LocalAssemblerInterface::getIntPtLiquidPressure);
    add_secondary_variable("mole_fraction_gas", 1,
                           &LocalAssemblerInterface::getIntPtMoleFractionGas);
    add_secondary_variable(
        "mole_fraction_liquid", 1,
        &LocalAssemblerInterface::getIntPtMoleFractionLiquid);
    add_secondary_variable("mass_fraction_gas", 1,
                           &LocalAssemblerInterface::getIntPtMassFractionGas);
    add_secondary_variable(
        "mass_fraction_liquid", 1,
        &LocalAssemblerInterface::getIntPtMassFractionLiquid);

    add_secondary_variable(
        "relative_permeability_gas", 1,
        &LocalAssemblerInterface::getIntPtRelativePermeabilityGas);
    add_secondary_variable(
        "relative_permeability_liquid", 1,
        &LocalAssemblerInterface::getIntPtRelativePermeabilityLiquid);

    add_secondary_variable("internal_energy_gas", 1,
                           &LocalAssemblerInterface::getIntPtInternalEnergyGas);

    add_secondary_variable(
        "internal_energy_liquid", 1,
        &LocalAssemblerInterface::getIntPtInternalEnergyLiquid);

    add_secondary_variable("enthalpy_gas", 1,
                           &LocalAssemblerInterface::getIntPtEnthalpyGas);

    add_secondary_variable("enthalpy_liquid", 1,
                           &LocalAssemblerInterface::getIntPtEnthalpyLiquid);

    add_secondary_variable("enthalpy_solid", 1,
                           &LocalAssemblerInterface::getIntPtEnthalpySolid);

    add_secondary_variable("enthalpy_CG", 1,
                           &LocalAssemblerInterface::getIntPtEnthalpyCG);

    add_secondary_variable("enthalpy_WG", 1,
                           &LocalAssemblerInterface::getIntPtEnthalpyWG);

    _process_data.element_saturation = MeshLib::getOrCreateMeshProperty<double>(
        const_cast<MeshLib::Mesh&>(mesh), "saturation_avg",
        MeshLib::MeshItemType::Cell, 1);

    _process_data.gas_pressure_interpolated =
        MeshLib::getOrCreateMeshProperty<double>(
            const_cast<MeshLib::Mesh&>(mesh), "gas_pressure_interpolated",
            MeshLib::MeshItemType::Node, 1);

    _process_data.capillary_pressure_interpolated =
        MeshLib::getOrCreateMeshProperty<double>(
            const_cast<MeshLib::Mesh&>(mesh), "capillary_pressure_interpolated",
            MeshLib::MeshItemType::Node, 1);

    _process_data.temperature_interpolated =
        MeshLib::getOrCreateMeshProperty<double>(
            const_cast<MeshLib::Mesh&>(mesh), "temperature_interpolated",
            MeshLib::MeshItemType::Node, 1);

    // Initialize local assemblers after all variables have been set.
    GlobalExecutor::executeMemberOnDereferenced(
        &LocalAssemblerInterface::initialize, _local_assemblers,
        *_local_to_global_index_map);
}

template <int DisplacementDim>
void TH2MProcess<DisplacementDim>::initializeBoundaryConditions()
{
    if (_use_monolithic_scheme)
    {
        initializeProcessBoundaryConditionsAndSourceTerms(
            *_local_to_global_index_map, monolithic_process_id);
        return;
    }

    // Staggered scheme:
    OGS_FATAL("A Staggered version of TH2M is not implemented.");
}

template <int DisplacementDim>
void TH2MProcess<DisplacementDim>::setInitialConditionsConcreteProcess(
    std::vector<GlobalVector*>& x, double const t, int const process_id)
{
    if (process_id != 0)
    {
        return;
    }

    DBUG("Set initial conditions of TH2MProcess.");

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];
    GlobalExecutor::executeSelectedMemberOnDereferenced(
        &LocalAssemblerInterface::setInitialConditions, _local_assemblers,
        pv.getActiveElementIDs(), getDOFTable(process_id), *x[process_id], t,
        _use_monolithic_scheme, process_id);
}

template <int DisplacementDim>
void TH2MProcess<DisplacementDim>::assembleConcreteProcess(
    const double t, double const dt, std::vector<GlobalVector*> const& x,
    std::vector<GlobalVector*> const& xdot, int const process_id,
    GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b)
{
    DBUG("Assemble the equations for TH2M");

    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_table = {std::ref(*_local_to_global_index_map)};
    // Call global assembler for each local assembly item.
    GlobalExecutor::executeMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assemble, _local_assemblers,
        dof_table, t, dt, x, xdot, process_id, M, K, b);
}

template <int DisplacementDim>
void TH2MProcess<DisplacementDim>::assembleWithJacobianConcreteProcess(
    const double t, double const dt, std::vector<GlobalVector*> const& x,
    std::vector<GlobalVector*> const& xdot, const double dxdot_dx,
    const double dx_dx, int const process_id, GlobalMatrix& M, GlobalMatrix& K,
    GlobalVector& b, GlobalMatrix& Jac)
{
    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_tables;
    // For the monolithic scheme
    if (_use_monolithic_scheme)
    {
        DBUG("Assemble the Jacobian of TH2M for the monolithic scheme.");
        dof_tables.emplace_back(*_local_to_global_index_map);
    }
    else
    {
        // For the staggered scheme
        OGS_FATAL("A Staggered version of TH2M is not implemented.");
    }

    GlobalExecutor::executeMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assembleWithJacobian,
        _local_assemblers, dof_tables, t, dt, x, xdot, dxdot_dx, dx_dx,
        process_id, M, K, b, Jac);

    auto copyRhs = [&](int const variable_id, auto& output_vector) {
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
void TH2MProcess<DisplacementDim>::preTimestepConcreteProcess(
    std::vector<GlobalVector*> const& x, double const t, double const dt,
    const int process_id)
{
    DBUG("PreTimestep TH2MProcess.");

    if (hasMechanicalProcess(process_id))
    {
        GlobalExecutor::executeMemberOnDereferenced(
            &LocalAssemblerInterface::preTimestep, _local_assemblers,
            *_local_to_global_index_map, *x[process_id], t, dt);
    }
}

template <int DisplacementDim>
void TH2MProcess<DisplacementDim>::postTimestepConcreteProcess(
    std::vector<GlobalVector*> const& x, double const t, double const dt,
    const int /*process_id*/)
{
    DBUG("PostTimestep TH2MProcess.");
    auto const dof_tables = getDOFTables(x.size());
    GlobalExecutor::executeMemberOnDereferenced(
        &LocalAssemblerInterface::postTimestep, _local_assemblers, dof_tables,
        x, t, dt);
}

template <int DisplacementDim>
void TH2MProcess<DisplacementDim>::postNonLinearSolverConcreteProcess(
    GlobalVector const& x, GlobalVector const& xdot, const double t,
    double const dt, const int process_id)
{
    if (!hasMechanicalProcess(process_id))
    {
        return;
    }

    DBUG("PostNonLinearSolver TH2MProcess.");
    // Calculate strain, stress or other internal variables of mechanics.
    GlobalExecutor::executeMemberOnDereferenced(
        &LocalAssemblerInterface::postNonLinearSolver, _local_assemblers,
        getDOFTable(process_id), x, xdot, t, dt, _use_monolithic_scheme,
        process_id);
}

template <int DisplacementDim>
void TH2MProcess<DisplacementDim>::computeSecondaryVariableConcrete(
    double const t, double const dt, std::vector<GlobalVector*> const& x,
    GlobalVector const& x_dot, const int process_id)
{
    if (process_id != 0)
    {
        return;
    }

    DBUG("Compute the secondary variables for TH2MProcess.");
    auto const dof_tables = getDOFTables(x.size());

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];
    GlobalExecutor::executeSelectedMemberOnDereferenced(
        &LocalAssemblerInterface::computeSecondaryVariable, _local_assemblers,
        pv.getActiveElementIDs(), dof_tables, t, dt, x, x_dot, process_id);
}

template <int DisplacementDim>
std::tuple<NumLib::LocalToGlobalIndexMap*, bool>
TH2MProcess<DisplacementDim>::getDOFTableForExtrapolatorData() const
{
    const bool manage_storage = false;
    return std::make_tuple(_local_to_global_index_map_single_component.get(),
                           manage_storage);
}

template <int DisplacementDim>
NumLib::LocalToGlobalIndexMap const& TH2MProcess<DisplacementDim>::getDOFTable(
    const int process_id) const
{
    if (hasMechanicalProcess(process_id))
    {
        return *_local_to_global_index_map;
    }

    // For the equation of pressure
    return *_local_to_global_index_map_with_base_nodes;
}

template <int DisplacementDim>
std::vector<NumLib::LocalToGlobalIndexMap const*>
TH2MProcess<DisplacementDim>::getDOFTables(int const number_of_processes) const
{
    std::vector<NumLib::LocalToGlobalIndexMap const*> dof_tables;
    dof_tables.reserve(number_of_processes);
    std::generate_n(std::back_inserter(dof_tables), number_of_processes,
                    [&]() { return &getDOFTable(dof_tables.size()); });
    return dof_tables;
}
template class TH2MProcess<2>;
template class TH2MProcess<3>;

}  // namespace TH2M
}  // namespace ProcessLib
