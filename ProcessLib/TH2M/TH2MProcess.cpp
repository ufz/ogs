/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "TH2MProcess.h"

#include <cassert>

#include "CreateTH2MLocalAssemblers.h"
#include "MaterialLib/SolidModels/MechanicsBase.h"  // for the instantiation of process data
#include "MathLib/KelvinVector.h"
#include "MeshLib/Elements/Utils.h"
#include "MeshLib/Utils/getOrCreateMeshProperty.h"
#include "NumLib/DOF/DOFTableUtil.h"
#include "ProcessLib/Deformation/SolidMaterialInternalToSecondaryVariables.h"
#include "ProcessLib/Process.h"
#include "ProcessLib/Reflection/ReflectionForExtrapolation.h"
#include "ProcessLib/Reflection/ReflectionForIPWriters.h"
#include "ProcessLib/Utils/ComputeResiduum.h"
#include "ProcessLib/Utils/SetIPDataInitialConditions.h"
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
      AssemblyMixin<TH2MProcess<DisplacementDim>>{*_jacobian_assembler},
      _process_data(std::move(process_data))
{
    ProcessLib::Reflection::addReflectedIntegrationPointWriters<
        DisplacementDim>(
        LocalAssemblerInterface<DisplacementDim>::getReflectionDataForOutput(),
        _integration_point_writer, integration_order, local_assemblers_);
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
    _mesh_subset_all_nodes = std::make_unique<MeshLib::MeshSubset>(
        _mesh, _mesh.getNodes(), _process_data.use_TaylorHood_elements);
    // Create single component dof in the mesh's base nodes.
    _base_nodes = MeshLib::getBaseNodes(_mesh.getElements());
    _mesh_subset_base_nodes = std::make_unique<MeshLib::MeshSubset>(
        _mesh, _base_nodes, _process_data.use_TaylorHood_elements);

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
    createLocalAssemblers<DisplacementDim>(
        mesh.getElements(), dof_table, local_assemblers_,
        NumLib::IntegrationOrder{integration_order}, mesh.isAxiallySymmetric(),
        _process_data);

    auto add_secondary_variable = [&](std::string const& name,
                                      int const num_components,
                                      auto get_ip_values_function)
    {
        _secondary_variables.addSecondaryVariable(
            name,
            makeExtrapolator(num_components, getExtrapolator(),
                             local_assemblers_,
                             std::move(get_ip_values_function)));
    };

    ProcessLib::Reflection::addReflectedSecondaryVariables<DisplacementDim>(
        LocalAssemblerInterface<DisplacementDim>::getReflectionDataForOutput(),
        _secondary_variables, getExtrapolator(), local_assemblers_);

    add_secondary_variable(
        "velocity_gas", mesh.getDimension(),
        &LocalAssemblerInterface<DisplacementDim>::getIntPtDarcyVelocityGas);
    add_secondary_variable(
        "velocity_liquid", mesh.getDimension(),
        &LocalAssemblerInterface<DisplacementDim>::getIntPtDarcyVelocityLiquid);
    add_secondary_variable(
        "diffusion_velocity_vapour_gas", mesh.getDimension(),
        &LocalAssemblerInterface<
            DisplacementDim>::getIntPtDiffusionVelocityVapourGas);
    add_secondary_variable(
        "diffusion_velocity_gas_gas", mesh.getDimension(),
        &LocalAssemblerInterface<
            DisplacementDim>::getIntPtDiffusionVelocityGasGas);
    add_secondary_variable(
        "diffusion_velocity_solute_liquid", mesh.getDimension(),
        &LocalAssemblerInterface<
            DisplacementDim>::getIntPtDiffusionVelocitySoluteLiquid);
    add_secondary_variable(
        "diffusion_velocity_liquid_liquid", mesh.getDimension(),
        &LocalAssemblerInterface<
            DisplacementDim>::getIntPtDiffusionVelocityLiquidLiquid);

    add_secondary_variable(
        "vapour_pressure", 1,
        &LocalAssemblerInterface<DisplacementDim>::getIntPtVapourPressure);
    add_secondary_variable(
        "porosity", 1,
        &LocalAssemblerInterface<DisplacementDim>::getIntPtPorosity);
    add_secondary_variable(
        "gas_density", 1,
        &LocalAssemblerInterface<DisplacementDim>::getIntPtGasDensity);
    add_secondary_variable(
        "solid_density", 1,
        &LocalAssemblerInterface<DisplacementDim>::getIntPtSolidDensity);
    add_secondary_variable(
        "liquid_density", 1,
        &LocalAssemblerInterface<DisplacementDim>::getIntPtLiquidDensity);
    add_secondary_variable(
        "mole_fraction_gas", 1,
        &LocalAssemblerInterface<DisplacementDim>::getIntPtMoleFractionGas);
    add_secondary_variable(
        "mass_fraction_gas", 1,
        &LocalAssemblerInterface<DisplacementDim>::getIntPtMassFractionGas);
    add_secondary_variable(
        "mass_fraction_liquid", 1,
        &LocalAssemblerInterface<DisplacementDim>::getIntPtMassFractionLiquid);

    add_secondary_variable(
        "enthalpy_gas", 1,
        &LocalAssemblerInterface<DisplacementDim>::getIntPtEnthalpyGas);

    add_secondary_variable(
        "enthalpy_liquid", 1,
        &LocalAssemblerInterface<DisplacementDim>::getIntPtEnthalpyLiquid);

    add_secondary_variable(
        "enthalpy_solid", 1,
        &LocalAssemblerInterface<DisplacementDim>::getIntPtEnthalpySolid);

    //
    // enable output of internal variables defined by material models
    //
    ProcessLib::Deformation::solidMaterialInternalToSecondaryVariables<
        LocalAssemblerInterface<DisplacementDim>>(_process_data.solid_materials,
                                                  add_secondary_variable);

    ProcessLib::Deformation::
        solidMaterialInternalVariablesToIntegrationPointWriter(
            _process_data.solid_materials, local_assemblers_,
            _integration_point_writer, integration_order);

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

    _process_data.liquid_pressure_interpolated =
        MeshLib::getOrCreateMeshProperty<double>(
            const_cast<MeshLib::Mesh&>(mesh), "liquid_pressure_interpolated",
            MeshLib::MeshItemType::Node, 1);

    _process_data.temperature_interpolated =
        MeshLib::getOrCreateMeshProperty<double>(
            const_cast<MeshLib::Mesh&>(mesh), "temperature_interpolated",
            MeshLib::MeshItemType::Node, 1);

    setIPDataInitialConditions(_integration_point_writer, mesh.getProperties(),
                               local_assemblers_);

    // Initialize local assemblers after all variables have been set.
    GlobalExecutor::executeMemberOnDereferenced(
        &LocalAssemblerInterface<DisplacementDim>::initialize,
        local_assemblers_, *_local_to_global_index_map);
}

template <int DisplacementDim>
void TH2MProcess<DisplacementDim>::initializeBoundaryConditions(
    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const& media)
{
    if (_use_monolithic_scheme)
    {
        initializeProcessBoundaryConditionsAndSourceTerms(
            *_local_to_global_index_map, monolithic_process_id, media);
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

    auto get_a_dof_table_func = [this](const int process_id) -> auto&
    {
        return getDOFTable(process_id);
    };
    GlobalExecutor::executeMemberOnDereferenced(
        &LocalAssemblerInterface<DisplacementDim>::setInitialConditions,
        local_assemblers_, NumLib::getDOFTables(x.size(), get_a_dof_table_func),
        x, t, process_id);
}

template <int DisplacementDim>
void TH2MProcess<DisplacementDim>::assembleConcreteProcess(
    const double t, double const dt, std::vector<GlobalVector*> const& x,
    std::vector<GlobalVector*> const& x_prev, int const process_id,
    GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b)
{
    DBUG("Assemble the equations for TH2M");

    AssemblyMixin<TH2MProcess<DisplacementDim>>::assemble(t, dt, x, x_prev,
                                                          process_id, M, K, b);
}

template <int DisplacementDim>
void TH2MProcess<DisplacementDim>::assembleWithJacobianConcreteProcess(
    const double t, double const dt, std::vector<GlobalVector*> const& x,
    std::vector<GlobalVector*> const& x_prev, int const process_id,
    GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b, GlobalMatrix& Jac)
{
    if (!_use_monolithic_scheme)
    {
        OGS_FATAL("A Staggered version of TH2M is not implemented.");
    }

    AssemblyMixin<TH2MProcess<DisplacementDim>>::assembleWithJacobian(
        t, dt, x, x_prev, process_id, M, K, b, Jac);
}

template <int DisplacementDim>
void TH2MProcess<DisplacementDim>::preTimestepConcreteProcess(
    std::vector<GlobalVector*> const& x, double const t, double const dt,
    const int process_id)
{
    DBUG("PreTimestep TH2MProcess.");

    if (hasMechanicalProcess(process_id))
    {
        ProcessLib::ProcessVariable const& pv =
            getProcessVariables(process_id)[0];

        GlobalExecutor::executeSelectedMemberOnDereferenced(
            &LocalAssemblerInterface<DisplacementDim>::preTimestep,
            local_assemblers_, pv.getActiveElementIDs(),
            *_local_to_global_index_map, *x[process_id], t, dt);
    }

    AssemblyMixin<TH2MProcess<DisplacementDim>>::updateActiveElements(
        process_id);
}

template <int DisplacementDim>
void TH2MProcess<DisplacementDim>::postTimestepConcreteProcess(
    std::vector<GlobalVector*> const& x,
    std::vector<GlobalVector*> const& x_prev, double const t, double const dt,
    const int process_id)
{
    DBUG("PostTimestep TH2MProcess.");

    auto get_a_dof_table_func = [this](const int processe_id) -> auto&
    {
        return getDOFTable(processe_id);
    };
    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];
    GlobalExecutor::executeSelectedMemberOnDereferenced(
        &LocalAssemblerInterface<DisplacementDim>::postTimestep,
        local_assemblers_, pv.getActiveElementIDs(),
        NumLib::getDOFTables(x.size(), get_a_dof_table_func), x, x_prev, t, dt,
        process_id);
}

template <int DisplacementDim>
void TH2MProcess<DisplacementDim>::computeSecondaryVariableConcrete(
    double const t, double const dt, std::vector<GlobalVector*> const& x,
    GlobalVector const& x_prev, const int process_id)
{
    if (process_id != 0)
    {
        return;
    }

    DBUG("Compute the secondary variables for TH2MProcess.");

    auto get_a_dof_table_func = [this](const int processe_id) -> auto&
    {
        return getDOFTable(processe_id);
    };
    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];
    GlobalExecutor::executeSelectedMemberOnDereferenced(
        &LocalAssemblerInterface<DisplacementDim>::computeSecondaryVariable,
        local_assemblers_, pv.getActiveElementIDs(),
        NumLib::getDOFTables(x.size(), get_a_dof_table_func), t, dt, x, x_prev,
        process_id);
}

template <int DisplacementDim>
std::vector<std::string>
TH2MProcess<DisplacementDim>::initializeAssemblyOnSubmeshes(
    std::vector<std::reference_wrapper<MeshLib::Mesh>> const& meshes)
{
    INFO("TH2M process initializeSubmeshOutput().");
    const int process_id = 0;
    std::vector<std::string> residuum_names{
        "GasMassFlowRate", "LiquidMassFlowRate", "HeatFlowRate", "NodalForces"};

    AssemblyMixin<TH2MProcess<DisplacementDim>>::initializeAssemblyOnSubmeshes(
        process_id, meshes, residuum_names);

    return residuum_names;
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

template class TH2MProcess<2>;
template class TH2MProcess<3>;

}  // namespace TH2M
}  // namespace ProcessLib
