/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
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
#include "ProcessLib/Deformation/SolidMaterialInternalToSecondaryVariables.h"
#include "ProcessLib/Reflection/ReflectionForExtrapolation.h"
#include "ProcessLib/Reflection/ReflectionForIPWriters.h"
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
      AssemblyMixin<RichardsMechanicsProcess<DisplacementDim>>{
          *_jacobian_assembler},
      process_data_(std::move(process_data))
{
    nodal_forces_ = MeshLib::getOrCreateMeshProperty<double>(
        mesh, "NodalForces", MeshLib::MeshItemType::Node, DisplacementDim);

    hydraulic_flow_ = MeshLib::getOrCreateMeshProperty<double>(
        mesh, "MassFlowRate", MeshLib::MeshItemType::Node, 1);

    ProcessLib::Reflection::addReflectedIntegrationPointWriters<
        DisplacementDim>(LocalAssemblerIF::getReflectionDataForOutput(),
                         _integration_point_writer, integration_order,
                         local_assemblers_);
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
    auto const& l = *local_to_global_index_map_with_base_nodes_;
    return {l.dofSizeWithoutGhosts(), l.dofSizeWithoutGhosts(),
            &l.getGhostIndices(), &sparsity_pattern_with_linear_element_};
}

template <int DisplacementDim>
void RichardsMechanicsProcess<DisplacementDim>::constructDofTable()
{
    // Create single component dof in every of the mesh's nodes.
    _mesh_subset_all_nodes =
        std::make_unique<MeshLib::MeshSubset>(_mesh, _mesh.getNodes());
    // Create single component dof in the mesh's base nodes.
    base_nodes_ = MeshLib::getBaseNodes(_mesh.getElements());
    mesh_subset_base_nodes_ =
        std::make_unique<MeshLib::MeshSubset>(_mesh, base_nodes_);

    // TODO move the two data members somewhere else.
    // for extrapolation of secondary variables of stress or strain
    std::vector<MeshLib::MeshSubset> all_mesh_subsets_single_component{
        *_mesh_subset_all_nodes};
    local_to_global_index_map_single_component_ =
        std::make_unique<NumLib::LocalToGlobalIndexMap>(
            std::move(all_mesh_subsets_single_component),
            // by location order is needed for output
            NumLib::ComponentOrder::BY_LOCATION);

    if (_use_monolithic_scheme)
    {
        // For pressure, which is the first
        std::vector<MeshLib::MeshSubset> all_mesh_subsets{
            *mesh_subset_base_nodes_};

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
            *mesh_subset_base_nodes_};
        local_to_global_index_map_with_base_nodes_ =
            std::make_unique<NumLib::LocalToGlobalIndexMap>(
                std::move(all_mesh_subsets_base_nodes),
                // by location order is needed for output
                NumLib::ComponentOrder::BY_LOCATION);

        sparsity_pattern_with_linear_element_ = NumLib::computeSparsityPattern(
            *local_to_global_index_map_with_base_nodes_, _mesh);

        assert(_local_to_global_index_map);
        assert(local_to_global_index_map_with_base_nodes_);
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
        mesh.getElements(), dof_table, local_assemblers_,
        NumLib::IntegrationOrder{integration_order}, mesh.isAxiallySymmetric(),
        process_data_);

    ProcessLib::Reflection::addReflectedSecondaryVariables<DisplacementDim>(
        LocalAssemblerIF::getReflectionDataForOutput(), _secondary_variables,
        getExtrapolator(), local_assemblers_);

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

    //
    // enable output of internal variables defined by material models
    //
    ProcessLib::Deformation::solidMaterialInternalToSecondaryVariables<
        LocalAssemblerIF>(process_data_.solid_materials,
                          add_secondary_variable);

    ProcessLib::Deformation::
        solidMaterialInternalVariablesToIntegrationPointWriter(
            process_data_.solid_materials, local_assemblers_,
            _integration_point_writer, integration_order);

    process_data_.element_saturation = MeshLib::getOrCreateMeshProperty<double>(
        const_cast<MeshLib::Mesh&>(mesh), "saturation_avg",
        MeshLib::MeshItemType::Cell, 1);

    process_data_.element_porosity = MeshLib::getOrCreateMeshProperty<double>(
        const_cast<MeshLib::Mesh&>(mesh), "porosity_avg",
        MeshLib::MeshItemType::Cell, 1);

    process_data_.element_stresses = MeshLib::getOrCreateMeshProperty<double>(
        const_cast<MeshLib::Mesh&>(mesh), "stress_avg",
        MeshLib::MeshItemType::Cell,
        MathLib::KelvinVector::KelvinVectorType<
            DisplacementDim>::RowsAtCompileTime);

    process_data_.pressure_interpolated =
        MeshLib::getOrCreateMeshProperty<double>(
            const_cast<MeshLib::Mesh&>(mesh), "pressure_interpolated",
            MeshLib::MeshItemType::Node, 1);

    setIPDataInitialConditions(_integration_point_writer, mesh.getProperties(),
                               local_assemblers_);

    // Initialize local assemblers after all variables have been set.
    GlobalExecutor::executeMemberOnDereferenced(&LocalAssemblerIF::initialize,
                                                local_assemblers_,
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
        *local_to_global_index_map_with_base_nodes_, hydraulic_process_id,
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

    GlobalExecutor::executeSelectedMemberOnDereferenced(
        &LocalAssemblerIF::setInitialConditions, local_assemblers_,
        getActiveElementIDs(), getDOFTables(x.size()), x, t, process_id);
}

template <int DisplacementDim>
void RichardsMechanicsProcess<DisplacementDim>::assembleConcreteProcess(
    const double t, double const dt, std::vector<GlobalVector*> const& x,
    std::vector<GlobalVector*> const& x_prev, int const process_id,
    GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b)
{
    DBUG("Assemble the equations for RichardsMechanics");

    AssemblyMixin<RichardsMechanicsProcess<DisplacementDim>>::assemble(
        t, dt, x, x_prev, process_id, M, K, b);
}

template <int DisplacementDim>
void RichardsMechanicsProcess<DisplacementDim>::
    assembleWithJacobianConcreteProcess(
        const double t, double const dt, std::vector<GlobalVector*> const& x,
        std::vector<GlobalVector*> const& x_prev, int const process_id,
        GlobalVector& b, GlobalMatrix& Jac)
{
    // For the monolithic scheme
    if (_use_monolithic_scheme)
    {
        DBUG(
            "Assemble the Jacobian of RichardsMechanics for the monolithic"
            " scheme.");
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
    }

    auto const dof_tables = getDOFTables(x.size());
    AssemblyMixin<RichardsMechanicsProcess<DisplacementDim>>::
        assembleWithJacobian(t, dt, x, x_prev, process_id, b, Jac);

    auto copyRhs = [&](int const variable_id, auto& output_vector)
    {
        if (_use_monolithic_scheme)
        {
            transformVariableFromGlobalVector(b, variable_id, *dof_tables[0],
                                              output_vector,
                                              std::negate<double>());
        }
        else
        {
            transformVariableFromGlobalVector(b, 0, *dof_tables[process_id],
                                              output_vector,
                                              std::negate<double>());
        }
    };
    if (_use_monolithic_scheme || process_id == 0)
    {
        copyRhs(0, *hydraulic_flow_);
    }
    if (_use_monolithic_scheme || process_id == 1)
    {
        copyRhs(1, *nodal_forces_);
    }
}

template <int DisplacementDim>
void RichardsMechanicsProcess<DisplacementDim>::preTimestepConcreteProcess(
    std::vector<GlobalVector*> const& x, double const t, double const dt,
    const int process_id)
{
    DBUG("PreTimestep RichardsMechanicsProcess.");

    GlobalExecutor::executeSelectedMemberOnDereferenced(
        &LocalAssemblerIF::preTimestep, local_assemblers_,
        getActiveElementIDs(), *_local_to_global_index_map, *x[process_id], t,
        dt);

    AssemblyMixin<
        RichardsMechanicsProcess<DisplacementDim>>::updateActiveElements();
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

        GlobalExecutor::executeSelectedMemberOnDereferenced(
            &LocalAssemblerIF::postTimestep, local_assemblers_,
            getActiveElementIDs(), getDOFTables(x.size()), x, x_prev, t, dt,
            process_id);
    }
}

template <int DisplacementDim>
std::vector<std::vector<std::string>>
RichardsMechanicsProcess<DisplacementDim>::initializeAssemblyOnSubmeshes(
    std::vector<std::reference_wrapper<MeshLib::Mesh>> const& meshes)
{
    INFO("RichardsMechanics process initializeSubmeshOutput().");
    std::vector<std::vector<std::string>> residuum_names{
        {"MassFlowRate", "NodalForces"}};

    AssemblyMixin<RichardsMechanicsProcess<DisplacementDim>>::
        initializeAssemblyOnSubmeshes(meshes, residuum_names);

    return residuum_names;
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

    GlobalExecutor::executeSelectedMemberOnDereferenced(
        &LocalAssemblerIF::computeSecondaryVariable, local_assemblers_,
        getActiveElementIDs(), getDOFTables(x.size()), t, dt, x, x_prev,
        process_id);
}

template <int DisplacementDim>
std::tuple<NumLib::LocalToGlobalIndexMap*, bool> RichardsMechanicsProcess<
    DisplacementDim>::getDOFTableForExtrapolatorData() const
{
    const bool manage_storage = false;
    return std::make_tuple(local_to_global_index_map_single_component_.get(),
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
    return *local_to_global_index_map_with_base_nodes_;
}

template class RichardsMechanicsProcess<2>;
template class RichardsMechanicsProcess<3>;

}  // namespace RichardsMechanics
}  // namespace ProcessLib
