/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ThermoRichardsMechanicsProcess.h"

#include <cassert>

#include "BaseLib/Error.h"
#include "MeshLib/Elements/Utils.h"
#include "NumLib/DOF/DOFTableUtil.h"
#include "ProcessLib/Deformation/SolidMaterialInternalToSecondaryVariables.h"
#include "ProcessLib/Process.h"
#include "ProcessLib/Reflection/ReflectionForExtrapolation.h"
#include "ProcessLib/Reflection/ReflectionForIPWriters.h"
#include "ProcessLib/Utils/CreateLocalAssemblersTaylorHood.h"
#include "ProcessLib/Utils/SetIPDataInitialConditions.h"
#include "ThermoRichardsMechanicsFEM.h"

namespace ProcessLib
{
namespace ThermoRichardsMechanics
{
template <typename ConstitutiveTraits>
struct ThermoRichardsMechanicsLocalAssembler3Args
{
    template <typename ShapeFunctionDisplacement, typename ShapeFunction,
              int DisplacementDim>
    using LocalAssemblerImplementation =
        ThermoRichardsMechanicsLocalAssembler<ShapeFunctionDisplacement,
                                              ShapeFunction, DisplacementDim,
                                              ConstitutiveTraits>;
};

template <int DisplacementDim, typename ConstitutiveTraits>
ThermoRichardsMechanicsProcess<DisplacementDim, ConstitutiveTraits>::
    ThermoRichardsMechanicsProcess(
        std::string name,
        MeshLib::Mesh& mesh,
        std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&&
            jacobian_assembler,
        std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const&
            parameters,
        unsigned const integration_order,
        std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
            process_variables,
        ThermoRichardsMechanicsProcessData<DisplacementDim,
                                           ConstitutiveTraits>&& process_data,
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
        mesh, "MassFlowRate", MeshLib::MeshItemType::Node, 1);

    heat_flux_ = MeshLib::getOrCreateMeshProperty<double>(
        mesh, "HeatFlowRate", MeshLib::MeshItemType::Node, 1);

    ProcessLib::Reflection::addReflectedIntegrationPointWriters<
        DisplacementDim>(LocalAssemblerIF::getReflectionDataForOutput(),
                         _integration_point_writer, integration_order,
                         local_assemblers_);
}

template <int DisplacementDim, typename ConstitutiveTraits>
bool ThermoRichardsMechanicsProcess<DisplacementDim,
                                    ConstitutiveTraits>::isLinear() const
{
    return false;
}

template <int DisplacementDim, typename ConstitutiveTraits>
MathLib::MatrixSpecifications ThermoRichardsMechanicsProcess<
    DisplacementDim,
    ConstitutiveTraits>::getMatrixSpecifications(const int /*process_id*/) const
{
    auto const& l = *_local_to_global_index_map;
    return {l.dofSizeWithoutGhosts(), l.dofSizeWithoutGhosts(),
            &l.getGhostIndices(), &this->_sparsity_pattern};
}

template <int DisplacementDim, typename ConstitutiveTraits>
void ThermoRichardsMechanicsProcess<DisplacementDim,
                                    ConstitutiveTraits>::constructDofTable()
{
    // Create single component dof in every of the mesh's nodes.
    _mesh_subset_all_nodes = std::make_unique<MeshLib::MeshSubset>(
        _mesh, _mesh.getNodes(), process_data_.use_TaylorHood_elements);
    // Create single component dof in the mesh's base nodes.
    base_nodes_ = MeshLib::getBaseNodes(_mesh.getElements());
    mesh_subset_base_nodes_ = std::make_unique<MeshLib::MeshSubset>(
        _mesh, base_nodes_, process_data_.use_TaylorHood_elements);

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

template <int DisplacementDim, typename ConstitutiveTraits>
void ThermoRichardsMechanicsProcess<DisplacementDim, ConstitutiveTraits>::
    initializeConcreteProcess(NumLib::LocalToGlobalIndexMap const& dof_table,
                              MeshLib::Mesh const& mesh,
                              unsigned const integration_order)
{
    createLocalAssemblersHM<
        DisplacementDim,
        ThermoRichardsMechanicsLocalAssembler3Args<
            ConstitutiveTraits>::template LocalAssemblerImplementation>(
        mesh.getElements(), dof_table, local_assemblers_,
        NumLib::IntegrationOrder{integration_order}, mesh.isAxiallySymmetric(),
        process_data_);

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
        LocalAssemblerIF::getReflectionDataForOutput(), _secondary_variables,
        getExtrapolator(), local_assemblers_);

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

    process_data_.element_liquid_density =
        MeshLib::getOrCreateMeshProperty<double>(
            const_cast<MeshLib::Mesh&>(mesh), "liquid_density_avg",
            MeshLib::MeshItemType::Cell, 1);

    process_data_.element_viscosity = MeshLib::getOrCreateMeshProperty<double>(
        const_cast<MeshLib::Mesh&>(mesh), "viscosity_avg",
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

    bool const remove_name_suffix = true;
    setIPDataInitialConditions(_integration_point_writer, mesh.getProperties(),
                               local_assemblers_, remove_name_suffix);

    // Initialize local assemblers after all variables have been set.
    GlobalExecutor::executeMemberOnDereferenced(&LocalAssemblerIF::initialize,
                                                local_assemblers_,
                                                *_local_to_global_index_map);
}

template <int DisplacementDim, typename ConstitutiveTraits>
void ThermoRichardsMechanicsProcess<
    DisplacementDim, ConstitutiveTraits>::initializeBoundaryConditions()
{
    const int process_id = 0;
    initializeProcessBoundaryConditionsAndSourceTerms(
        *_local_to_global_index_map, process_id);
}

template <int DisplacementDim, typename ConstitutiveTraits>
void ThermoRichardsMechanicsProcess<DisplacementDim, ConstitutiveTraits>::
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

template <int DisplacementDim, typename ConstitutiveTraits>
void ThermoRichardsMechanicsProcess<DisplacementDim, ConstitutiveTraits>::
    assembleConcreteProcess(const double /*t*/, double const /*dt*/,
                            std::vector<GlobalVector*> const& /*x*/,
                            std::vector<GlobalVector*> const& /*xdot*/,
                            int const /*process_id*/, GlobalMatrix& /*M*/,
                            GlobalMatrix& /*K*/, GlobalVector& /*b*/)
{
    OGS_FATAL(
        "The Picard method or the Newton-Raphson method with numerical "
        "Jacobian is not implemented for ThermoRichardsMechanics with the full "
        "monolithic coupling scheme");
}

template <int DisplacementDim, typename ConstitutiveTraits>
void ThermoRichardsMechanicsProcess<DisplacementDim, ConstitutiveTraits>::
    assembleWithJacobianConcreteProcess(const double t, double const dt,
                                        std::vector<GlobalVector*> const& x,
                                        std::vector<GlobalVector*> const& xdot,
                                        int const process_id, GlobalMatrix& M,
                                        GlobalMatrix& K, GlobalVector& b,
                                        GlobalMatrix& Jac)
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
        process_id, M, K, b, Jac);

    auto copyRhs = [&](int const variable_id, auto& output_vector)
    {
        transformVariableFromGlobalVector(b, variable_id, dof_tables[0],
                                          output_vector, std::negate<double>());
    };

    copyRhs(0, *heat_flux_);
    copyRhs(1, *hydraulic_flow_);
    copyRhs(2, *nodal_forces_);
}

template <int DisplacementDim, typename ConstitutiveTraits>
void ThermoRichardsMechanicsProcess<DisplacementDim, ConstitutiveTraits>::
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

template <int DisplacementDim, typename ConstitutiveTraits>
void ThermoRichardsMechanicsProcess<DisplacementDim, ConstitutiveTraits>::
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

template <int DisplacementDim, typename ConstitutiveTraits>
std::tuple<NumLib::LocalToGlobalIndexMap*, bool> ThermoRichardsMechanicsProcess<
    DisplacementDim, ConstitutiveTraits>::getDOFTableForExtrapolatorData() const
{
    const bool manage_storage = false;
    return std::make_tuple(local_to_global_index_map_single_component_.get(),
                           manage_storage);
}

template <int DisplacementDim, typename ConstitutiveTraits>
NumLib::LocalToGlobalIndexMap const& ThermoRichardsMechanicsProcess<
    DisplacementDim, ConstitutiveTraits>::getDOFTable(const int /*process_id*/)
    const
{
    return *_local_to_global_index_map;
}

template <int DisplacementDim, typename ConstitutiveTraits>
std::vector<NumLib::LocalToGlobalIndexMap const*>
ThermoRichardsMechanicsProcess<DisplacementDim, ConstitutiveTraits>::
    getDOFTables(int const number_of_processes) const
{
    std::vector<NumLib::LocalToGlobalIndexMap const*> dof_tables;
    dof_tables.reserve(number_of_processes);
    std::generate_n(std::back_inserter(dof_tables), number_of_processes,
                    [&]() { return &getDOFTable(dof_tables.size()); });
    return dof_tables;
}

template class ThermoRichardsMechanicsProcess<
    2, ConstitutiveStress_StrainTemperature::ConstitutiveTraits<2>>;
template class ThermoRichardsMechanicsProcess<
    3, ConstitutiveStress_StrainTemperature::ConstitutiveTraits<3>>;

#if OGS_USE_MFRONT
template class ThermoRichardsMechanicsProcess<
    2, ConstitutiveStressSaturation_StrainPressureTemperature::
           ConstitutiveTraits<2>>;
template class ThermoRichardsMechanicsProcess<
    3, ConstitutiveStressSaturation_StrainPressureTemperature::
           ConstitutiveTraits<3>>;
#endif

}  // namespace ThermoRichardsMechanics
}  // namespace ProcessLib
