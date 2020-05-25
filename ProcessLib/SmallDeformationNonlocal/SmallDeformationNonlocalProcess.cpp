/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "SmallDeformationNonlocalProcess.h"

#include <nlohmann/json.hpp>
#include <iostream>

#include "ProcessLib/Output/IntegrationPointWriter.h"

// Reusing local assembler creation code.
#include "ProcessLib/SmallDeformation/CreateLocalAssemblers.h"

namespace ProcessLib
{
namespace SmallDeformationNonlocal
{
template <int DisplacementDim>
SmallDeformationNonlocalProcess<DisplacementDim>::
    SmallDeformationNonlocalProcess(
        std::string name,
        MeshLib::Mesh& mesh,
        std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&&
            jacobian_assembler,
        std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const&
            parameters,
        unsigned const integration_order,
        std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
            process_variables,
        SmallDeformationNonlocalProcessData<DisplacementDim>&& process_data,
        SecondaryVariableCollection&& secondary_variables)
    : Process(std::move(name), mesh, std::move(jacobian_assembler), parameters,
              integration_order, std::move(process_variables),
              std::move(secondary_variables)),
      process_data_(std::move(process_data))
{
    nodal_forces_ = MeshLib::getOrCreateMeshProperty<double>(
        mesh, "NodalForces", MeshLib::MeshItemType::Node, DisplacementDim);

    integration_point_writer_.emplace_back(
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

    integration_point_writer_.emplace_back(
        std::make_unique<IntegrationPointWriter>(
            "kappa_d_ip", 1 /*n_components*/, integration_order, [this]() {
                // Result containing integration point data for each local
                // assembler.
                std::vector<std::vector<double>> result;
                result.resize(local_assemblers_.size());

                for (std::size_t i = 0; i < local_assemblers_.size(); ++i)
                {
                    auto const& local_asm = *local_assemblers_[i];

                    result[i] = local_asm.getKappaD();
                }

                return result;
            }));
}

template <int DisplacementDim>
bool SmallDeformationNonlocalProcess<DisplacementDim>::isLinear() const
{
    return false;
}

template <int DisplacementDim>
void SmallDeformationNonlocalProcess<DisplacementDim>::
    initializeConcreteProcess(NumLib::LocalToGlobalIndexMap const& dof_table,
                              MeshLib::Mesh const& mesh,
                              unsigned const integration_order)
{
    // Reusing local assembler creation code.
    ProcessLib::SmallDeformation::createLocalAssemblers<
        DisplacementDim, SmallDeformationNonlocalLocalAssembler>(
        mesh.getElements(), dof_table, local_assemblers_,
        mesh.isAxiallySymmetric(), integration_order, process_data_);

    // TODO move the two data members somewhere else.
    // for extrapolation of secondary variables
    std::vector<MeshLib::MeshSubset> all_mesh_subsets_single_component{
        *mesh_subset_all_nodes_};
    local_to_global_index_map_single_component_ =
        std::make_unique<NumLib::LocalToGlobalIndexMap>(
            std::move(all_mesh_subsets_single_component),
            // by location order is needed for output
            NumLib::ComponentOrder::BY_LOCATION);

    Process::secondary_variables_.addSecondaryVariable(
        "sigma",
        makeExtrapolator(MathLib::KelvinVector::KelvinVectorType<
                             DisplacementDim>::RowsAtCompileTime,
                         getExtrapolator(), local_assemblers_,
                         &LocalAssemblerInterface::getIntPtSigma));

    Process::secondary_variables_.addSecondaryVariable(
        "epsilon",
        makeExtrapolator(MathLib::KelvinVector::KelvinVectorType<
                             DisplacementDim>::RowsAtCompileTime,
                         getExtrapolator(), local_assemblers_,
                         &LocalAssemblerInterface::getIntPtEpsilon));

    Process::secondary_variables_.addSecondaryVariable(
        "eps_p_V",
        makeExtrapolator(1, getExtrapolator(), local_assemblers_,
                         &LocalAssemblerInterface::getIntPtEpsPV));
    Process::secondary_variables_.addSecondaryVariable(
        "eps_p_D_xx",
        makeExtrapolator(1, getExtrapolator(), local_assemblers_,
                         &LocalAssemblerInterface::getIntPtEpsPDXX));

    Process::secondary_variables_.addSecondaryVariable(
        "damage",
        makeExtrapolator(1, getExtrapolator(), local_assemblers_,
                         &LocalAssemblerInterface::getIntPtDamage));

    GlobalExecutor::executeMemberOnDereferenced(
        &LocalAssemblerInterface::nonlocal, local_assemblers_,
        local_assemblers_);

    // Set initial conditions for integration point data.
    for (auto const& ip_writer : integration_point_writer_)
    {
        auto const& name = ip_writer->name();
        // First check the field data, which is used for restart.
        if (mesh.getProperties().existsPropertyVector<double>(name))
        {
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
                mesh_property.getNumberOfComponents())
            {
                OGS_FATAL(
                    "Different number of components in meta data ({:d}) than "
                    "in "
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
        else if (mesh.getProperties().existsPropertyVector<double>(name +
                                                                   "_ic"))
        {  // Try to find cell data with '_ic' suffix
            auto const& mesh_property =
                *mesh.getProperties().template getPropertyVector<double>(name +
                                                                         "_ic");
            if (mesh_property.getMeshItemType() != MeshLib::MeshItemType::Cell)
            {
                continue;
            }

            // Now we have a vtk's cell data array containing the initial
            // conditions for the corresponding integration point writer.

            // For each assembler use the single cell value for all
            // integration points.
            for (std::size_t i = 0; i < local_assemblers_.size(); ++i)
            {
                auto& local_asm = local_assemblers_[i];

                std::vector<double> value(
                    &mesh_property[i],
                    &mesh_property[i] + mesh_property.getNumberOfComponents());
                // TODO (naumov) Check sizes / read size / etc.
                // OR reconstruct dimensions from size / component =
                // ip_points
                local_asm->setIPDataInitialConditionsFromCellData(name, value);
            }
        }
    }

    // Initialize local assemblers after all variables have been set.
    GlobalExecutor::executeMemberOnDereferenced(
        &LocalAssemblerInterface::initialize, local_assemblers_,
        *local_to_global_index_map_);
}

template <int DisplacementDim>
void SmallDeformationNonlocalProcess<DisplacementDim>::assembleConcreteProcess(
    const double t, double const dt, std::vector<GlobalVector*> const& x,
    std::vector<GlobalVector*> const& xdot, int const process_id,
    GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b)
{
    DBUG("Assemble SmallDeformationNonlocalProcess.");

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
void SmallDeformationNonlocalProcess<
    DisplacementDim>::preAssembleConcreteProcess(const double t,
                                                 double const dt,
                                                 GlobalVector const& x)
{
    DBUG("preAssemble SmallDeformationNonlocalProcess.");

    const int process_id = 0;
    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    // Call global assembler for each local assembly item.
    GlobalExecutor::executeSelectedMemberDereferenced(
        global_assembler_, &VectorMatrixAssembler::preAssemble,
        local_assemblers_, pv.getActiveElementIDs(),
        *local_to_global_index_map_, t, dt, x);
}

template <int DisplacementDim>
void SmallDeformationNonlocalProcess<DisplacementDim>::
    assembleWithJacobianConcreteProcess(
        const double t, double const dt, std::vector<GlobalVector*> const& x,
        GlobalVector const& xdot, const double dxdot_dx, const double dx_dx,
        int const process_id, GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b,
        GlobalMatrix& Jac)
{
    DBUG("AssembleWithJacobian SmallDeformationNonlocalProcess.");

    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_table = {std::ref(*local_to_global_index_map_)};

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    // Call global assembler for each local assembly item.
    GlobalExecutor::executeSelectedMemberDereferenced(
        global_assembler_, &VectorMatrixAssembler::assembleWithJacobian,
        local_assemblers_, pv.getActiveElementIDs(), dof_table, t, dt, x, xdot,
        dxdot_dx, dx_dx, process_id, M, K, b, Jac, coupled_solutions_);

    b.copyValues(*nodal_forces_);
    std::transform(nodal_forces_->begin(), nodal_forces_->end(),
                   nodal_forces_->begin(), [](double val) { return -val; });
}

template <int DisplacementDim>
void SmallDeformationNonlocalProcess<DisplacementDim>::
    postTimestepConcreteProcess(std::vector<GlobalVector*> const& x,
                                double const t,
                                double const dt,
                                int const process_id)
{
    DBUG("PostTimestep SmallDeformationNonlocalProcess.");

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    GlobalExecutor::executeSelectedMemberOnDereferenced(
        &LocalAssemblerInterface::postTimestep, local_assemblers_,
        pv.getActiveElementIDs(), *local_to_global_index_map_, *x[process_id],
        t, dt);
}

template <int DisplacementDim>
NumLib::IterationResult
SmallDeformationNonlocalProcess<DisplacementDim>::postIterationConcreteProcess(
    GlobalVector const& x)
{
    process_data_.crack_volume_old = process_data_.crack_volume;
    process_data_.crack_volume = 0.0;

    DBUG("PostNonLinearSolver crack volume computation.");

    const int process_id = 0;
    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    GlobalExecutor::executeSelectedMemberOnDereferenced(
        &LocalAssemblerInterface::computeCrackIntegral, local_assemblers_,
        pv.getActiveElementIDs(), *local_to_global_index_map_, x,
        process_data_.crack_volume);

    INFO("Integral of crack: {:g}", process_data_.crack_volume);

    return NumLib::IterationResult::SUCCESS;
}

template class SmallDeformationNonlocalProcess<2>;
template class SmallDeformationNonlocalProcess<3>;

}  // namespace SmallDeformationNonlocal
}  // namespace ProcessLib
