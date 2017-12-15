/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "SmallDeformationNonlocalProcess.h"

#include <nlohmann/json.hpp>
#include <iostream>

// Reusing local assembler creation code.
#include "ProcessLib/SmallDeformation/CreateLocalAssemblers.h"

namespace ProcessLib
{
namespace SmallDeformationNonlocal
{
template <int DisplacementDim>
SmallDeformationNonlocalProcess<DisplacementDim>::
    SmallDeformationNonlocalProcess(
        MeshLib::Mesh& mesh,
        std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&&
            jacobian_assembler,
        std::vector<std::unique_ptr<ParameterBase>> const& parameters,
        unsigned const integration_order,
        std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
            process_variables,
        SmallDeformationNonlocalProcessData<DisplacementDim>&& process_data,
        SecondaryVariableCollection&& secondary_variables,
        NumLib::NamedFunctionCaller&& named_function_caller)
    : Process(mesh, std::move(jacobian_assembler), parameters,
              integration_order, std::move(process_variables),
              std::move(secondary_variables), std::move(named_function_caller)),
      _process_data(std::move(process_data))
{
    _nodal_forces = MeshLib::getOrCreateMeshProperty<double>(
        mesh, "NodalForces", MeshLib::MeshItemType::Node, DisplacementDim);

    _integration_point_writer.emplace_back(
        std::make_unique<SigmaIntegrationPointWriter>(
            static_cast<int>(mesh.getDimension() == 2 ? 4 : 6) /*n components*/,
            integration_order, [this]() {
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

    _integration_point_writer.emplace_back(
        std::make_unique<KappaDIntegrationPointWriter>(
            integration_order, [this]() {
                // Result containing integration point data for each local
                // assembler.
                std::vector<std::vector<double>> result;
                result.resize(_local_assemblers.size());

                for (std::size_t i = 0; i < _local_assemblers.size(); ++i)
                {
                    auto const& local_asm = *_local_assemblers[i];

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
        mesh.getElements(), dof_table, _local_assemblers,
        mesh.isAxiallySymmetric(), integration_order, _process_data);

    // TODO move the two data members somewhere else.
    // for extrapolation of secondary variables
    std::vector<MeshLib::MeshSubset> all_mesh_subsets_single_component{
        *_mesh_subset_all_nodes};
    _local_to_global_index_map_single_component =
        std::make_unique<NumLib::LocalToGlobalIndexMap>(
            std::move(all_mesh_subsets_single_component),
            // by location order is needed for output
            NumLib::ComponentOrder::BY_LOCATION);

    Process::_secondary_variables.addSecondaryVariable(
        "sigma",
        makeExtrapolator(MathLib::KelvinVector::KelvinVectorType<
                             DisplacementDim>::RowsAtCompileTime,
                         getExtrapolator(), _local_assemblers,
                         &LocalAssemblerInterface::getIntPtSigma));

    Process::_secondary_variables.addSecondaryVariable(
        "epsilon",
        makeExtrapolator(MathLib::KelvinVector::KelvinVectorType<
                             DisplacementDim>::RowsAtCompileTime,
                         getExtrapolator(), _local_assemblers,
                         &LocalAssemblerInterface::getIntPtEpsilon));

    Process::_secondary_variables.addSecondaryVariable(
        "eps_p_V",
        makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                         &LocalAssemblerInterface::getIntPtEpsPV));
    Process::_secondary_variables.addSecondaryVariable(
        "eps_p_D_xx",
        makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                         &LocalAssemblerInterface::getIntPtEpsPDXX));

    Process::_secondary_variables.addSecondaryVariable(
        "damage",
        makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                         &LocalAssemblerInterface::getIntPtDamage));

    GlobalExecutor::executeMemberOnDereferenced(
        &LocalAssemblerInterface::nonlocal, _local_assemblers,
        _local_assemblers);


    // Set initial conditions for integration point data.
    for (auto const& ip_writer : _integration_point_writer)
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
                    "Different number of components in meta data (%d) than in "
                    "the integration point field data for '%s': %d.",
                    ip_meta_data.n_components, name.c_str(),
                    mesh_property.getNumberOfComponents());
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
            for (std::size_t i = 0; i < _local_assemblers.size(); ++i)
            {
                auto& local_asm = _local_assemblers[i];

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
}

template <int DisplacementDim>
void SmallDeformationNonlocalProcess<DisplacementDim>::assembleConcreteProcess(
    const double t, GlobalVector const& x, GlobalMatrix& M, GlobalMatrix& K,
    GlobalVector& b)
{
    DBUG("Assemble SmallDeformationNonlocalProcess.");

    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_table = {std::ref(*_local_to_global_index_map)};
    // Call global assembler for each local assembly item.
    GlobalExecutor::executeMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assemble, _local_assemblers,
        dof_table, t, x, M, K, b, _coupled_solutions);
}

template <int DisplacementDim>
void SmallDeformationNonlocalProcess<
    DisplacementDim>::preAssembleConcreteProcess(const double t,
                                                 GlobalVector const& x)
{
    DBUG("preAssemble SmallDeformationNonlocalProcess.");

    // Call global assembler for each local assembly item.
    GlobalExecutor::executeMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::preAssemble,
        _local_assemblers, *_local_to_global_index_map, t, x);
}

template <int DisplacementDim>
void SmallDeformationNonlocalProcess<DisplacementDim>::
    assembleWithJacobianConcreteProcess(const double t, GlobalVector const& x,
                                        GlobalVector const& xdot,
                                        const double dxdot_dx,
                                        const double dx_dx, GlobalMatrix& M,
                                        GlobalMatrix& K, GlobalVector& b,
                                        GlobalMatrix& Jac)
{
    DBUG("AssembleWithJacobian SmallDeformationNonlocalProcess.");

    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_table = {std::ref(*_local_to_global_index_map)};
    // Call global assembler for each local assembly item.
    GlobalExecutor::executeMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assembleWithJacobian,
        _local_assemblers, dof_table, t, x, xdot, dxdot_dx, dx_dx, M, K, b, Jac,
        _coupled_solutions);

    b.copyValues(*_nodal_forces);
    std::transform(_nodal_forces->begin(), _nodal_forces->end(),
                   _nodal_forces->begin(), [](double val) { return -val; });
}

template <int DisplacementDim>
void SmallDeformationNonlocalProcess<
    DisplacementDim>::preTimestepConcreteProcess(GlobalVector const& x,
                                                 double const t,
                                                 double const dt,
                                                 int const /*process_id*/)
{
    DBUG("PreTimestep SmallDeformationNonlocalProcess.");

    _process_data.dt = dt;
    _process_data.t = t;

    GlobalExecutor::executeMemberOnDereferenced(
        &LocalAssemblerInterface::preTimestep, _local_assemblers,
        *_local_to_global_index_map, x, t, dt);
}

template <int DisplacementDim>
NumLib::IterationResult
SmallDeformationNonlocalProcess<DisplacementDim>::postIterationConcreteProcess(
    GlobalVector const& x)
{
    _process_data.crack_volume_old = _process_data.crack_volume;
    _process_data.crack_volume = 0.0;

    DBUG("PostNonLinearSolver crack volume computation.");

    GlobalExecutor::executeMemberOnDereferenced(
        &LocalAssemblerInterface::computeCrackIntegral, _local_assemblers,
        *_local_to_global_index_map, x, _process_data.crack_volume);

    INFO("Integral of crack: %g", _process_data.crack_volume);

    return NumLib::IterationResult::SUCCESS;
}

template class SmallDeformationNonlocalProcess<2>;
template class SmallDeformationNonlocalProcess<3>;

}  // namespace SmallDeformationNonlocal
}  // namespace ProcessLib
