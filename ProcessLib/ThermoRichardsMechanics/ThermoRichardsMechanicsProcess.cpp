/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ThermoRichardsMechanicsProcess.h"

#include <cassert>

#include "BaseLib/Error.h"
#include "CreateLocalAssemblers.h"
#include "MeshLib/Elements/Utils.h"
#include "NumLib/DOF/ComputeSparsityPattern.h"
#include "ProcessLib/Deformation/SolidMaterialInternalToSecondaryVariables.h"
#include "ProcessLib/Process.h"
#include "ThermoRichardsMechanicsFEM.h"

namespace ProcessLib
{
namespace ThermoRichardsMechanics
{
template <int DisplacementDim>
ThermoRichardsMechanicsProcess<DisplacementDim>::ThermoRichardsMechanicsProcess(
    std::string name,
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    unsigned const integration_order,
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
        process_variables,
    ThermoRichardsMechanicsProcessData<DisplacementDim>&& process_data,
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

    // TODO (naumov) remove ip suffix. Probably needs modification of the mesh
    // properties, s.t. there is no "overlapping" with cell/point data.
    // See getOrCreateMeshProperty.
    _integration_point_writer.emplace_back(
        std::make_unique<IntegrationPointWriter>(
            "sigma_ip",
            static_cast<int>(mesh.getDimension() == 2 ? 4 : 6) /*n components*/,
            integration_order, [this]() {
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

    _integration_point_writer.emplace_back(
        std::make_unique<IntegrationPointWriter>(
            "saturation_ip", 1 /*n components*/, integration_order, [this]() {
                // Result containing integration point data for each local
                // assembler.
                std::vector<std::vector<double>> result;
                result.resize(local_assemblers_.size());

                for (std::size_t i = 0; i < local_assemblers_.size(); ++i)
                {
                    auto const& local_asm = *local_assemblers_[i];
                    result[i] = local_asm.getSaturation();
                }

                return result;
            }));

    _integration_point_writer.emplace_back(
        std::make_unique<IntegrationPointWriter>(
            "porosity_ip", 1 /*n components*/, integration_order, [this]() {
                // Result containing integration point data for each local
                // assembler.
                std::vector<std::vector<double>> result;
                result.resize(local_assemblers_.size());

                for (std::size_t i = 0; i < local_assemblers_.size(); ++i)
                {
                    auto const& local_asm = *local_assemblers_[i];
                    result[i] = local_asm.getPorosity();
                }

                return result;
            }));

    _integration_point_writer.emplace_back(
        std::make_unique<IntegrationPointWriter>(
            "transport_porosity_ip", 1 /*n components*/, integration_order,
            [this]() {
                // Result containing integration point data for each local
                // assembler.
                std::vector<std::vector<double>> result;
                result.resize(local_assemblers_.size());

                for (std::size_t i = 0; i < local_assemblers_.size(); ++i)
                {
                    auto const& local_asm = *local_assemblers_[i];
                    result[i] = local_asm.getTransportPorosity();
                }

                return result;
            }));

    _integration_point_writer.emplace_back(
        std::make_unique<IntegrationPointWriter>(
            "swelling_stress_ip",
            static_cast<int>(mesh.getDimension() == 2 ? 4 : 6) /*n components*/,
            integration_order, [this]() {
                // Result containing integration point data for each local
                // assembler.
                std::vector<std::vector<double>> result;
                result.resize(local_assemblers_.size());

                for (std::size_t i = 0; i < local_assemblers_.size(); ++i)
                {
                    auto const& local_asm = *local_assemblers_[i];
                    result[i] = local_asm.getSwellingStress();
                }

                return result;
            }));

    _integration_point_writer.emplace_back(
        std::make_unique<IntegrationPointWriter>(
            "epsilon_ip",
            static_cast<int>(mesh.getDimension() == 2 ? 4 : 6) /*n components*/,
            integration_order, [this]() {
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
bool ThermoRichardsMechanicsProcess<DisplacementDim>::isLinear() const
{
    return false;
}

template <int DisplacementDim>
MathLib::MatrixSpecifications
ThermoRichardsMechanicsProcess<DisplacementDim>::getMatrixSpecifications(
    const int /*process_id*/) const
{
    auto const& l = *_local_to_global_index_map;
    return {l.dofSizeWithoutGhosts(), l.dofSizeWithoutGhosts(),
            &l.getGhostIndices(), &this->_sparsity_pattern};
}

template <int DisplacementDim>
void ThermoRichardsMechanicsProcess<DisplacementDim>::constructDofTable()
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
    std::vector<MeshLib::MeshSubset> all__meshsubsets_single_component{
        *_mesh_subset_all_nodes};
    local_to_global_index_map_single_component_ =
        std::make_unique<NumLib::LocalToGlobalIndexMap>(
            std::move(all__meshsubsets_single_component),
            // by location order is needed for output
            NumLib::ComponentOrder::BY_LOCATION);

    // For temperature, which is the first variable.
    std::vector<MeshLib::MeshSubset> all__meshsubsets{*mesh_subset_base_nodes_};

    // For pressure, which is the second variable
    all__meshsubsets.push_back(*mesh_subset_base_nodes_);

    // For displacement.
    const int monolithic_process_id = 0;
    std::generate_n(std::back_inserter(all__meshsubsets),
                    getProcessVariables(monolithic_process_id)[2]
                        .get()
                        .getNumberOfGlobalComponents(),
                    [&]() { return *_mesh_subset_all_nodes; });

    std::vector<int> const vec_n_components{1, 1, DisplacementDim};
    _local_to_global_index_map =
        std::make_unique<NumLib::LocalToGlobalIndexMap>(
            std::move(all__meshsubsets), vec_n_components,
            NumLib::ComponentOrder::BY_LOCATION);
    assert(_local_to_global_index_map);
}

template <int DisplacementDim>
void ThermoRichardsMechanicsProcess<DisplacementDim>::initializeConcreteProcess(
    NumLib::LocalToGlobalIndexMap const& dof_table,
    MeshLib::Mesh const& mesh,
    unsigned const integration_order)
{
    using nlohmann::json;

    const int mechanical_process_id = 0;
    const int deformation_variable_id = 2;
    createLocalAssemblers<DisplacementDim,
                          ThermoRichardsMechanicsLocalAssembler>(
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
        _secondary_variables.addSecondaryVariable(
            name,
            makeExtrapolator(num_components, getExtrapolator(),
                             local_assemblers_,
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

    add_secondary_variable("porosity", 1, &LocalAssemblerIF::getIntPtPorosity);

    add_secondary_variable("transport_porosity", 1,
                           &LocalAssemblerIF::getIntPtTransportPorosity);

    add_secondary_variable("dry_density_solid", 1,
                           &LocalAssemblerIF::getIntPtDryDensitySolid);

    //
    // enable output of internal variables defined by material models
    //
    ProcessLib::Deformation::solidMaterialInternalToSecondaryVariables<
        LocalAssemblerIF>(process_data_.solid_materials,
                          add_secondary_variable);

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
    process_data_.temperature_interpolated =
        MeshLib::getOrCreateMeshProperty<double>(
            const_cast<MeshLib::Mesh&>(mesh), "temperature_interpolated",
            MeshLib::MeshItemType::Node, 1);

    // Set initial conditions for integration point data.
    for (auto const& ip_writer : _integration_point_writer)
    {
        // Find the mesh property with integration point writer's name.
        auto const& name = ip_writer->name();
        if (!mesh.getProperties().existsPropertyVector<double>(name))
        {
            continue;
        }
        auto const& _meshproperty =
            *mesh.getProperties().template getPropertyVector<double>(name);

        // The mesh property must be defined on integration points.
        if (_meshproperty.getMeshItemType() !=
            MeshLib::MeshItemType::IntegrationPoint)
        {
            continue;
        }

        auto const ip_meta_data = getIntegrationPointMetaData(mesh, name);

        // Check the number of components.
        if (ip_meta_data.n_components !=
            _meshproperty.getNumberOfGlobalComponents())
        {
            OGS_FATAL(
                "Different number of components in meta data ({:d}) than in "
                "the integration point field data for '{:s}': {:d}.",
                ip_meta_data.n_components, name,
                _meshproperty.getNumberOfGlobalComponents());
        }

        // Now we have a properly named vtk's field data array and the
        // corresponding meta data.
        std::size_t position = 0;
        for (auto& local_asm : local_assemblers_)
        {
            std::size_t const integration_points_read =
                local_asm->setIPDataInitialConditions(
                    name, &_meshproperty[position],
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
    GlobalExecutor::executeMemberOnDereferenced(&LocalAssemblerIF::initialize,
                                                local_assemblers_,
                                                *_local_to_global_index_map);
}

template <int DisplacementDim>
void ThermoRichardsMechanicsProcess<
    DisplacementDim>::initializeBoundaryConditions()
{
    const int process_id = 0;
    initializeProcessBoundaryConditionsAndSourceTerms(
        *_local_to_global_index_map, process_id);
}

template <int DisplacementDim>
void ThermoRichardsMechanicsProcess<DisplacementDim>::
    setInitialConditionsConcreteProcess(std::vector<GlobalVector*>& x,
                                        double const t,
                                        int const process_id)
{
    DBUG("SetInitialConditions ThermoRichardsMechanicsProcess.");

    GlobalExecutor::executeMemberOnDereferenced(
        &LocalAssemblerIF::setInitialConditions, local_assemblers_,
        *_local_to_global_index_map, *x[process_id], t, _use_monolithic_scheme,
        process_id);
}

template <int DisplacementDim>
void ThermoRichardsMechanicsProcess<DisplacementDim>::assembleConcreteProcess(
    const double /*t*/, double const /*dt*/,
    std::vector<GlobalVector*> const& /*x*/,
    std::vector<GlobalVector*> const& /*xdot*/, int const /*process_id*/,
    GlobalMatrix& /*M*/, GlobalMatrix& /*K*/, GlobalVector& /*b*/)
{
    OGS_FATAL(
        "The Picard method or the Newton-Raphson method with numerical "
        "Jacobian is not implemented for ThermoRichardsMechanics with the full "
        "monolithic coupling scheme");
}

template <int DisplacementDim>
void ThermoRichardsMechanicsProcess<DisplacementDim>::
    assembleWithJacobianConcreteProcess(
        const double t, double const dt, std::vector<GlobalVector*> const& x,
        std::vector<GlobalVector*> const& xdot, const double dxdot_dx,
        const double dx_dx, int const process_id, GlobalMatrix& M,
        GlobalMatrix& K, GlobalVector& b, GlobalMatrix& Jac)
{
    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_tables;

    DBUG(
        "Assemble the Jacobian of ThermoRichardsMechanics for the monolithic "
        "scheme.");
    dof_tables.emplace_back(*_local_to_global_index_map);

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    GlobalExecutor::executeSelectedMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assembleWithJacobian,
        local_assemblers_, pv.getActiveElementIDs(), dof_tables, t, dt, x, xdot,
        dxdot_dx, dx_dx, process_id, M, K, b, Jac);

    auto copyRhs = [&](int const variable_id, auto& output_vector) {
        transformVariableFromGlobalVector(b, variable_id, dof_tables[0],
                                          output_vector, std::negate<double>());
    };

    copyRhs(1, *hydraulic_flow_);
    copyRhs(2, *nodal_forces_);
}

template <int DisplacementDim>
void ThermoRichardsMechanicsProcess<DisplacementDim>::
    postTimestepConcreteProcess(std::vector<GlobalVector*> const& x,
                                double const t, double const dt,
                                const int process_id)
{
    DBUG("PostTimestep ThermoRichardsMechanicsProcess.");

    auto const dof_tables = getDOFTables(x.size());

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];
    GlobalExecutor::executeSelectedMemberOnDereferenced(
        &LocalAssemblerIF::postTimestep, local_assemblers_,
        pv.getActiveElementIDs(), dof_tables, x, t, dt);
}

template <int DisplacementDim>
void ThermoRichardsMechanicsProcess<DisplacementDim>::
    computeSecondaryVariableConcrete(const double t, const double dt,
                                     std::vector<GlobalVector*> const& x,
                                     GlobalVector const& x_dot,
                                     int const process_id)
{
    DBUG("Compute the secondary variables for ThermoRichardsMechanicsProcess.");

    auto const dof_tables = getDOFTables(x.size());

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    GlobalExecutor::executeSelectedMemberOnDereferenced(
        &LocalAssemblerIF::computeSecondaryVariable, local_assemblers_,
        pv.getActiveElementIDs(), dof_tables, t, dt, x, x_dot, process_id);
}

template <int DisplacementDim>
std::tuple<NumLib::LocalToGlobalIndexMap*, bool> ThermoRichardsMechanicsProcess<
    DisplacementDim>::getDOFTableForExtrapolatorData() const
{
    const bool manage_storage = false;
    return std::make_tuple(local_to_global_index_map_single_component_.get(),
                           manage_storage);
}

template <int DisplacementDim>
NumLib::LocalToGlobalIndexMap const&
ThermoRichardsMechanicsProcess<DisplacementDim>::getDOFTable(
    const int /*process_id*/) const
{
    return *_local_to_global_index_map;
}

template <int DisplacementDim>
std::vector<NumLib::LocalToGlobalIndexMap const*>
ThermoRichardsMechanicsProcess<DisplacementDim>::getDOFTables(
    int const number_of_processes) const
{
    std::vector<NumLib::LocalToGlobalIndexMap const*> dof_tables;
    dof_tables.reserve(number_of_processes);
    std::generate_n(std::back_inserter(dof_tables), number_of_processes,
                    [&]() { return &getDOFTable(dof_tables.size()); });
    return dof_tables;
}

template class ThermoRichardsMechanicsProcess<2>;
template class ThermoRichardsMechanicsProcess<3>;

}  // namespace ThermoRichardsMechanics
}  // namespace ProcessLib
