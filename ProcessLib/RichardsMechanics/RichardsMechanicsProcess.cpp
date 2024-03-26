/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "RichardsMechanicsProcess.h"

#include <cassert>

#include "MeshLib/Elements/Utils.h"
#include "MeshLib/Utils/getOrCreateMeshProperty.h"
#include "NumLib/DOF/ComputeSparsityPattern.h"
#include "NumLib/DOF/DOFTableUtil.h"
#include "ProcessLib/Deformation/SolidMaterialInternalToSecondaryVariables.h"
#include "ProcessLib/Utils/CreateLocalAssemblersTaylorHood.h"
#include "ProcessLib/Utils/SetIPDataInitialConditions.h"
#include "RichardsMechanicsFEM.h"
#include "RichardsMechanicsProcessData.h"

namespace ProcessLib
{
namespace RichardsMechanics
{
template <int DisplacementDim>
RichardsMechanicsProcess<DisplacementDim>::RichardsMechanicsProcess(
    std::string name,
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    unsigned const integration_order,
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
        process_variables,
    RichardsMechanicsProcessData<DisplacementDim>&& process_data,
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

    // TODO (naumov) remove ip suffix. Probably needs modification of the mesh
    // properties, s.t. there is no "overlapping" with cell/point data.
    // See getOrCreateMeshProperty.
    _integration_point_writer.emplace_back(
        std::make_unique<MeshLib::IntegrationPointWriter>(
            "sigma_ip",
            static_cast<int>(mesh.getDimension() == 2 ? 4 : 6) /*n components*/,
            integration_order, _local_assemblers, &LocalAssemblerIF::getSigma));

    _integration_point_writer.emplace_back(
        std::make_unique<MeshLib::IntegrationPointWriter>(
            "saturation_ip", 1 /*n components*/, integration_order,
            _local_assemblers, &LocalAssemblerIF::getSaturation));

    _integration_point_writer.emplace_back(
        std::make_unique<MeshLib::IntegrationPointWriter>(
            "porosity_ip", 1 /*n components*/, integration_order,
            _local_assemblers, &LocalAssemblerIF::getPorosity));

    _integration_point_writer.emplace_back(
        std::make_unique<MeshLib::IntegrationPointWriter>(
            "transport_porosity_ip", 1 /*n components*/, integration_order,
            _local_assemblers, &LocalAssemblerIF::getTransportPorosity));

    _integration_point_writer.emplace_back(
        std::make_unique<MeshLib::IntegrationPointWriter>(
            "swelling_stress_ip",
            static_cast<int>(mesh.getDimension() == 2 ? 4 : 6) /*n components*/,
            integration_order, _local_assemblers,
            &LocalAssemblerIF::getSwellingStress));

    _integration_point_writer.emplace_back(
        std::make_unique<MeshLib::IntegrationPointWriter>(
            "epsilon_ip",
            static_cast<int>(mesh.getDimension() == 2 ? 4 : 6) /*n components*/,
            integration_order, _local_assemblers,
            &LocalAssemblerIF::getEpsilon));
}

template <int DisplacementDim>
bool RichardsMechanicsProcess<DisplacementDim>::isLinear() const
{
    return false;
}

template <int DisplacementDim>
MathLib::MatrixSpecifications
RichardsMechanicsProcess<DisplacementDim>::getMatrixSpecifications(
    const int process_id) const
{
    // For the monolithic scheme or the M process (deformation) in the staggered
    // scheme.
    if (_use_monolithic_scheme || process_id == 1)
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
void RichardsMechanicsProcess<DisplacementDim>::constructDofTable()
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
void RichardsMechanicsProcess<DisplacementDim>::initializeConcreteProcess(
    NumLib::LocalToGlobalIndexMap const& dof_table,
    MeshLib::Mesh const& mesh,
    unsigned const integration_order)
{
    ProcessLib::createLocalAssemblersHM<DisplacementDim,
                                        RichardsMechanicsLocalAssembler>(
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

    add_secondary_variable("swelling_stress",
                           MathLib::KelvinVector::KelvinVectorType<
                               DisplacementDim>::RowsAtCompileTime,
                           &LocalAssemblerIF::getIntPtSwellingStress);

    add_secondary_variable("epsilon",
                           MathLib::KelvinVector::KelvinVectorType<
                               DisplacementDim>::RowsAtCompileTime,
                           &LocalAssemblerIF::getIntPtEpsilon);

    add_secondary_variable("velocity", DisplacementDim,
                           &LocalAssemblerIF::getIntPtDarcyVelocity);

    add_secondary_variable("saturation", 1,
                           &LocalAssemblerIF::getIntPtSaturation);

    add_secondary_variable("micro_saturation", 1,
                           &LocalAssemblerIF::getIntPtMicroSaturation);

    add_secondary_variable("micro_pressure", 1,
                           &LocalAssemblerIF::getIntPtMicroPressure);

    add_secondary_variable("porosity", 1, &LocalAssemblerIF::getIntPtPorosity);

    add_secondary_variable("transport_porosity", 1,
                           &LocalAssemblerIF::getIntPtTransportPorosity);

    add_secondary_variable("dry_density_solid", 1,
                           &LocalAssemblerIF::getIntPtDryDensitySolid);

    //
    // enable output of internal variables defined by material models
    //
    ProcessLib::Deformation::solidMaterialInternalToSecondaryVariables<
        LocalAssemblerIF>(_process_data.solid_materials,
                          add_secondary_variable);

    ProcessLib::Deformation::
        solidMaterialInternalVariablesToIntegrationPointWriter(
            _process_data.solid_materials, _local_assemblers,
            _integration_point_writer, integration_order);

    _process_data.element_saturation = MeshLib::getOrCreateMeshProperty<double>(
        const_cast<MeshLib::Mesh&>(mesh), "saturation_avg",
        MeshLib::MeshItemType::Cell, 1);

    _process_data.element_porosity = MeshLib::getOrCreateMeshProperty<double>(
        const_cast<MeshLib::Mesh&>(mesh), "porosity_avg",
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

    setIPDataInitialConditions(_integration_point_writer, mesh.getProperties(),
                               _local_assemblers);

    // Initialize local assemblers after all variables have been set.
    GlobalExecutor::executeMemberOnDereferenced(&LocalAssemblerIF::initialize,
                                                _local_assemblers,
                                                *_local_to_global_index_map);
}

template <int DisplacementDim>
void RichardsMechanicsProcess<DisplacementDim>::initializeBoundaryConditions(
    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const& media)
{
    if (_use_monolithic_scheme)
    {
        const int monolithic_process_id = 0;
        initializeProcessBoundaryConditionsAndSourceTerms(
            *_local_to_global_index_map, monolithic_process_id, media);
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
void RichardsMechanicsProcess<DisplacementDim>::
    setInitialConditionsConcreteProcess(std::vector<GlobalVector*>& x,
                                        double const t,
                                        int const process_id)
{
    if (process_id != 0)
    {
        return;
    }

    DBUG("SetInitialConditions RichardsMechanicsProcess.");

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    auto get_a_dof_table_func = [this](const int process_id) -> auto&
    {
        return getDOFTable(process_id);
    };
    GlobalExecutor::executeSelectedMemberOnDereferenced(
        &LocalAssemblerIF::setInitialConditions, _local_assemblers,
        pv.getActiveElementIDs(),
        NumLib::getDOFTables(x.size(), get_a_dof_table_func), x, t, process_id);
}

template <int DisplacementDim>
void RichardsMechanicsProcess<DisplacementDim>::assembleConcreteProcess(
    const double t, double const dt, std::vector<GlobalVector*> const& x,
    std::vector<GlobalVector*> const& x_prev, int const process_id,
    GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b)
{
    DBUG("Assemble the equations for RichardsMechanics");

    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_table = {std::ref(*_local_to_global_index_map)};
    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    // Call global assembler for each local assembly item.
    GlobalExecutor::executeSelectedMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assemble, _local_assemblers,
        pv.getActiveElementIDs(), dof_table, t, dt, x, x_prev, process_id, M, K,
        b);
}

template <int DisplacementDim>
void RichardsMechanicsProcess<DisplacementDim>::
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
            "Assemble the Jacobian of RichardsMechanics for the monolithic"
            " scheme.");
        dof_tables.emplace_back(*_local_to_global_index_map);
    }
    else
    {
        // For the staggered scheme
        if (process_id == 0)
        {
            DBUG(
                "Assemble the Jacobian equations of liquid fluid process in "
                "RichardsMechanics for the staggered scheme.");
        }
        else
        {
            DBUG(
                "Assemble the Jacobian equations of mechanical process in "
                "RichardsMechanics for the staggered scheme.");
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
        copyRhs(0, *_hydraulic_flow);
    }
    if (_use_monolithic_scheme || process_id == 1)
    {
        copyRhs(1, *_nodal_forces);
    }
}

template <int DisplacementDim>
void RichardsMechanicsProcess<DisplacementDim>::postTimestepConcreteProcess(
    std::vector<GlobalVector*> const& x,
    std::vector<GlobalVector*> const& x_prev, double const t, double const dt,
    const int process_id)
{
    if (hasMechanicalProcess(process_id))
    {
        DBUG("PostTimestep RichardsMechanicsProcess.");

        auto get_a_dof_table_func = [this](const int processe_id) -> auto&
        {
            return getDOFTable(processe_id);
        };
        ProcessLib::ProcessVariable const& pv =
            getProcessVariables(process_id)[0];
        GlobalExecutor::executeSelectedMemberOnDereferenced(
            &LocalAssemblerIF::postTimestep, _local_assemblers,
            pv.getActiveElementIDs(),
            NumLib::getDOFTables(x.size(), get_a_dof_table_func), x, x_prev, t,
            dt, process_id);
    }
}

template <int DisplacementDim>
void RichardsMechanicsProcess<DisplacementDim>::
    computeSecondaryVariableConcrete(const double t, const double dt,
                                     std::vector<GlobalVector*> const& x,
                                     GlobalVector const& x_prev,
                                     int const process_id)
{
    if (process_id != 0)
    {
        return;
    }

    DBUG("Compute the secondary variables for RichardsMechanicsProcess.");

    auto get_a_dof_table_func = [this](const int processe_id) -> auto&
    {
        return getDOFTable(processe_id);
    };
    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];
    GlobalExecutor::executeSelectedMemberOnDereferenced(
        &LocalAssemblerIF::computeSecondaryVariable, _local_assemblers,
        pv.getActiveElementIDs(),
        NumLib::getDOFTables(x.size(), get_a_dof_table_func), t, dt, x, x_prev,
        process_id);
}

template <int DisplacementDim>
std::tuple<NumLib::LocalToGlobalIndexMap*, bool> RichardsMechanicsProcess<
    DisplacementDim>::getDOFTableForExtrapolatorData() const
{
    const bool manage_storage = false;
    return std::make_tuple(_local_to_global_index_map_single_component.get(),
                           manage_storage);
}

template <int DisplacementDim>
NumLib::LocalToGlobalIndexMap const&
RichardsMechanicsProcess<DisplacementDim>::getDOFTable(
    const int process_id) const
{
    if (hasMechanicalProcess(process_id))
    {
        return *_local_to_global_index_map;
    }

    // For the equation of pressure
    return *_local_to_global_index_map_with_base_nodes;
}

template class RichardsMechanicsProcess<2>;
template class RichardsMechanicsProcess<3>;

}  // namespace RichardsMechanics
}  // namespace ProcessLib
