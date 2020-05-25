/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ThermoMechanicsProcess.h"

#include <cassert>

#include "NumLib/DOF/ComputeSparsityPattern.h"
#include "NumLib/DOF/DOFTableUtil.h"
#include "ProcessLib/Deformation/SolidMaterialInternalToSecondaryVariables.h"
#include "ProcessLib/Output/IntegrationPointWriter.h"
#include "ProcessLib/SmallDeformation/CreateLocalAssemblers.h"
#include "ThermoMechanicsFEM.h"

namespace ProcessLib
{
namespace ThermoMechanics
{
template <int DisplacementDim>
ThermoMechanicsProcess<DisplacementDim>::ThermoMechanicsProcess(
    std::string name, MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    unsigned const integration_order,
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
        process_variables,
    ThermoMechanicsProcessData<DisplacementDim>&& process_data,
    SecondaryVariableCollection&& secondary_variables,
    bool const use_monolithic_scheme)
    : Process(std::move(name), mesh, std::move(jacobian_assembler), parameters,
              integration_order, std::move(process_variables),
              std::move(secondary_variables), use_monolithic_scheme),
      process_data_(std::move(process_data))
{
    nodal_forces_ = MeshLib::getOrCreateMeshProperty<double>(
        mesh, "NodalForces", MeshLib::MeshItemType::Node, DisplacementDim);

    heat_flux_ = MeshLib::getOrCreateMeshProperty<double>(
        mesh, "HeatFlux", MeshLib::MeshItemType::Node, 1);

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

    integration_point_writer_.emplace_back(
        std::make_unique<IntegrationPointWriter>(
            "epsilon_m_ip",
            static_cast<int>(mesh.getDimension() == 2 ? 4 : 6) /*n components*/,
            2 /*integration order*/, [this]() {
                // Result containing integration point data for each local
                // assembler.
                std::vector<std::vector<double>> result;
                result.resize(local_assemblers_.size());

                for (std::size_t i = 0; i < local_assemblers_.size(); ++i)
                {
                    auto const& local_asm = *local_assemblers_[i];

                    result[i] = local_asm.getEpsilonMechanical();
                }

                return result;
            }));
}

template <int DisplacementDim>
bool ThermoMechanicsProcess<DisplacementDim>::isLinear() const
{
    return false;
}

template <int DisplacementDim>
MathLib::MatrixSpecifications
ThermoMechanicsProcess<DisplacementDim>::getMatrixSpecifications(
    const int process_id) const
{
    // For the monolithic scheme or the M process (deformation) in the staggered
    // scheme.
    if (use_monolithic_scheme_ ||
        process_id == process_data_.mechanics_process_id)
    {
        auto const& l = *local_to_global_index_map_;
        return {l.dofSizeWithoutGhosts(), l.dofSizeWithoutGhosts(),
                &l.getGhostIndices(), &this->sparsity_pattern_};
    }

    // For staggered scheme and T process.
    auto const& l = *local_to_global_index_map_single_component_;
    return {l.dofSizeWithoutGhosts(), l.dofSizeWithoutGhosts(),
            &l.getGhostIndices(), &sparsity_pattern_with_single_component_};
}

// TODO [WW]: remove if (use_monolithic_scheme_) during the refactoring of the
// coupling part.
template <int DisplacementDim>
void ThermoMechanicsProcess<DisplacementDim>::constructDofTable()
{
    // Note: the heat conduction process and the mechanical process use the same
    // order of shape functions.

    if (use_monolithic_scheme_)
    {
        constructMonolithicProcessDofTable();
        return;
    }
    constructDofTableOfSpecifiedProsessStaggerdScheme(
        process_data_.mechanics_process_id);

    // TODO move the two data members somewhere else.
    // for extrapolation of secondary variables of stress or strain
    std::vector<MeshLib::MeshSubset> all_mesh_subsets_single_component{
        *mesh_subset_all_nodes_};
    local_to_global_index_map_single_component_.reset(
        new NumLib::LocalToGlobalIndexMap(
            std::move(all_mesh_subsets_single_component),
            // by location order is needed for output
            NumLib::ComponentOrder::BY_LOCATION));

    if (!use_monolithic_scheme_)
    {
        sparsity_pattern_with_single_component_ =
            NumLib::computeSparsityPattern(
                *local_to_global_index_map_single_component_, mesh_);
    }
}

template <int DisplacementDim>
void ThermoMechanicsProcess<DisplacementDim>::initializeConcreteProcess(
    NumLib::LocalToGlobalIndexMap const& dof_table,
    MeshLib::Mesh const& mesh,
    unsigned const integration_order)
{
    ProcessLib::SmallDeformation::createLocalAssemblers<
        DisplacementDim, ThermoMechanicsLocalAssembler>(
        mesh.getElements(), dof_table, local_assemblers_,
        mesh.isAxiallySymmetric(), integration_order, process_data_);

    auto add_secondary_variable = [&](std::string const& name,
                                      int const num_components,
                                      auto get_ip_values_function) {
        secondary_variables_.addSecondaryVariable(
            name,
            makeExtrapolator(num_components, getExtrapolator(),
                             local_assemblers_,
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

    //
    // enable output of internal variables defined by material models
    //
    ProcessLib::Deformation::solidMaterialInternalToSecondaryVariables<
        LocalAssemblerInterface>(process_data_.solid_materials,
                                 add_secondary_variable);

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
        &LocalAssemblerInterface::initialize, local_assemblers_,
        *local_to_global_index_map_);
}

template <int DisplacementDim>
void ThermoMechanicsProcess<DisplacementDim>::initializeBoundaryConditions()
{
    if (use_monolithic_scheme_)
    {
        const int process_id_of_thermomechanics = 0;
        initializeProcessBoundaryConditionsAndSourceTerms(
            *local_to_global_index_map_, process_id_of_thermomechanics);
        return;
    }

    // Staggered scheme:
    // for the equations of heat conduction
    initializeProcessBoundaryConditionsAndSourceTerms(
        *local_to_global_index_map_single_component_,
        process_data_.heat_conduction_process_id);

    // for the equations of deformation.
    initializeProcessBoundaryConditionsAndSourceTerms(
        *local_to_global_index_map_, process_data_.mechanics_process_id);
}

template <int DisplacementDim>
void ThermoMechanicsProcess<DisplacementDim>::assembleConcreteProcess(
    const double t, double const dt, std::vector<GlobalVector*> const& x,
    std::vector<GlobalVector*> const& xdot, int const process_id,
    GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b)
{
    DBUG("Assemble ThermoMechanicsProcess.");

    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_table = {std::ref(*local_to_global_index_map_)};
    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    // Call global assembler for each local assembly item.
    GlobalExecutor::executeSelectedMemberDereferenced(
        global_assembler_, &VectorMatrixAssembler::assemble, local_assemblers_,
        pv.getActiveElementIDs(), dof_table, t, dt, x, xdot, process_id, M, K,
        b, coupled_solutions_);
}

template <int DisplacementDim>
void ThermoMechanicsProcess<DisplacementDim>::
    assembleWithJacobianConcreteProcess(
        const double t, double const dt, std::vector<GlobalVector*> const& x,
        GlobalVector const& xdot, const double dxdot_dx, const double dx_dx,
        int const process_id, GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b,
        GlobalMatrix& Jac)
{
    DBUG("AssembleJacobian ThermoMechanicsProcess.");

    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_tables;
    // For the monolithic scheme
    if (use_monolithic_scheme_)
    {
        DBUG(
            "Assemble the Jacobian of ThermoMechanics for the monolithic"
            " scheme.");
        dof_tables.emplace_back(*local_to_global_index_map_);
    }
    else
    {
        // For the staggered scheme
        if (process_id == process_data_.heat_conduction_process_id)
        {
            DBUG(
                "Assemble the Jacobian equations of heat conduction process in "
                "ThermoMechanics for the staggered scheme.");
        }
        else
        {
            DBUG(
                "Assemble the Jacobian equations of mechanical process in "
                "ThermoMechanics for the staggered scheme.");
        }

        // For the flexible appearance order of processes in the coupling.
        if (process_data_.heat_conduction_process_id ==
            0)  // First: the heat conduction process
        {
            dof_tables.emplace_back(
                *local_to_global_index_map_single_component_);
            dof_tables.emplace_back(*local_to_global_index_map_);
        }
        else  // vice versa
        {
            dof_tables.emplace_back(*local_to_global_index_map_);
            dof_tables.emplace_back(
                *local_to_global_index_map_single_component_);
        }

        setCoupledSolutionsOfPreviousTimeStep();
    }

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    GlobalExecutor::executeSelectedMemberDereferenced(
        global_assembler_, &VectorMatrixAssembler::assembleWithJacobian,
        local_assemblers_, pv.getActiveElementIDs(), dof_tables, t, dt, x, xdot,
        dxdot_dx, dx_dx, process_id, M, K, b, Jac, coupled_solutions_);

    // TODO (naumov): Refactor the copy rhs part. This is copy from HM.
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
    if (use_monolithic_scheme_ ||
        process_id == process_data_.heat_conduction_process_id)
    {
        copyRhs(0, *heat_flux_);
    }
    if (use_monolithic_scheme_ ||
        process_id == process_data_.mechanics_process_id)
    {
        copyRhs(1, *nodal_forces_);
    }
}

template <int DisplacementDim>
void ThermoMechanicsProcess<DisplacementDim>::preTimestepConcreteProcess(
    std::vector<GlobalVector*> const& x, double const t, double const dt,
    const int process_id)
{
    DBUG("PreTimestep ThermoMechanicsProcess.");

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    assert(process_id < 2);

    if (process_id == process_data_.mechanics_process_id)
    {
        GlobalExecutor::executeSelectedMemberOnDereferenced(
            &LocalAssemblerInterface::preTimestep, local_assemblers_,
            pv.getActiveElementIDs(), *local_to_global_index_map_,
            *x[process_id], t, dt);
        return;
    }

    // For the staggered scheme.
    if (!previous_T_)
    {
        previous_T_ = MathLib::MatrixVectorTraits<GlobalVector>::newInstance(
            *x[process_id]);
    }
    else
    {
        auto& x0 = *previous_T_;
        MathLib::LinAlg::copy(*x[process_id], x0);
    }

    auto& x0 = *previous_T_;
    MathLib::LinAlg::setLocalAccessibleVector(x0);
}

template <int DisplacementDim>
void ThermoMechanicsProcess<DisplacementDim>::postTimestepConcreteProcess(
    std::vector<GlobalVector*> const& x, double const t, double const dt,
    int const process_id)
{
    if (process_id != process_data_.mechanics_process_id)
    {
        return;
    }

    DBUG("PostTimestep ThermoMechanicsProcess.");

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    GlobalExecutor::executeSelectedMemberOnDereferenced(
        &LocalAssemblerInterface::postTimestep, local_assemblers_,
        pv.getActiveElementIDs(), *local_to_global_index_map_, *x[process_id],
        t, dt);
}

template <int DisplacementDim>
void ThermoMechanicsProcess<
    DisplacementDim>::setCoupledSolutionsOfPreviousTimeStep()
{
    coupled_solutions_->coupled_xs_t0.resize(1);
    coupled_solutions_->coupled_xs_t0[0] = previous_T_.get();
}

template <int DisplacementDim>
NumLib::LocalToGlobalIndexMap const&
ThermoMechanicsProcess<DisplacementDim>::getDOFTable(const int process_id) const
{
    if (process_data_.mechanics_process_id == process_id)
    {
        return *local_to_global_index_map_;
    }

    // For the equation of pressure
    return *local_to_global_index_map_single_component_;
}

template class ThermoMechanicsProcess<2>;
template class ThermoMechanicsProcess<3>;

}  // namespace ThermoMechanics
}  // namespace ProcessLib
