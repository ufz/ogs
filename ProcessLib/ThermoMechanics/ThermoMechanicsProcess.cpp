/**
 * \file
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ThermoMechanicsProcess.h"

#include <cassert>

#include "BaseLib/Functional.h"
#include "NumLib/DOF/ComputeSparsityPattern.h"
#include "NumLib/DOF/DOFTableUtil.h"
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
    NumLib::NamedFunctionCaller&& named_function_caller,
    bool const use_monolithic_scheme)
    : Process(std::move(name), mesh, std::move(jacobian_assembler), parameters,
              integration_order, std::move(process_variables),
              std::move(secondary_variables), std::move(named_function_caller),
              use_monolithic_scheme),
      _process_data(std::move(process_data))
{
    _nodal_forces = MeshLib::getOrCreateMeshProperty<double>(
        mesh, "NodalForces", MeshLib::MeshItemType::Node, DisplacementDim);

    _heat_flux = MeshLib::getOrCreateMeshProperty<double>(
        mesh, "HeatFlux", MeshLib::MeshItemType::Node, 1);

    _integration_point_writer.emplace_back(
        std::make_unique<IntegrationPointWriter>(
            "sigma_ip",
            static_cast<int>(mesh.getDimension() == 2 ? 4 : 6) /*n components*/,
            2 /*integration order*/, [this]() {
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
        std::make_unique<IntegrationPointWriter>(
            "epsilon_ip",
            static_cast<int>(mesh.getDimension() == 2 ? 4 : 6) /*n components*/,
            2 /*integration order*/, [this]() {
                // Result containing integration point data for each local
                // assembler.
                std::vector<std::vector<double>> result;
                result.resize(_local_assemblers.size());

                for (std::size_t i = 0; i < _local_assemblers.size(); ++i)
                {
                    auto const& local_asm = *_local_assemblers[i];

                    result[i] = local_asm.getEpsilon();
                }

                return result;
            }));

    _integration_point_writer.emplace_back(
        std::make_unique<IntegrationPointWriter>(
            "epsilon_m_ip",
            static_cast<int>(mesh.getDimension() == 2 ? 4 : 6) /*n components*/,
            2 /*integration order*/, [this]() {
                // Result containing integration point data for each local
                // assembler.
                std::vector<std::vector<double>> result;
                result.resize(_local_assemblers.size());

                for (std::size_t i = 0; i < _local_assemblers.size(); ++i)
                {
                    auto const& local_asm = *_local_assemblers[i];

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
    if (_use_monolithic_scheme ||
        process_id == _process_data.mechanics_process_id)
    {
        auto const& l = *_local_to_global_index_map;
        return {l.dofSizeWithoutGhosts(), l.dofSizeWithoutGhosts(),
                &l.getGhostIndices(), &this->_sparsity_pattern};
    }

    // For staggered scheme and T process.
    auto const& l = *_local_to_global_index_map_single_component;
    return {l.dofSizeWithoutGhosts(), l.dofSizeWithoutGhosts(),
            &l.getGhostIndices(), &_sparsity_pattern_with_single_component};
}

// TODO [WW]: remove if (_use_monolithic_scheme) during the refactoring of the
// coupling part.
template <int DisplacementDim>
void ThermoMechanicsProcess<DisplacementDim>::constructDofTable()
{
    // Note: the heat conduction process and the mechanical process use the same
    // order of shape functions.

    if (_use_monolithic_scheme)
    {
        constructMonolithicProcessDofTable();
        return;
    }
    constructDofTableOfSpecifiedProsessStaggerdScheme(
        _process_data.mechanics_process_id);

    // TODO move the two data members somewhere else.
    // for extrapolation of secondary variables of stress or strain
    std::vector<MeshLib::MeshSubset> all_mesh_subsets_single_component{
        *_mesh_subset_all_nodes};
    _local_to_global_index_map_single_component.reset(
        new NumLib::LocalToGlobalIndexMap(
            std::move(all_mesh_subsets_single_component),
            // by location order is needed for output
            NumLib::ComponentOrder::BY_LOCATION));

    if (!_use_monolithic_scheme)
    {
        _sparsity_pattern_with_single_component =
            NumLib::computeSparsityPattern(
                *_local_to_global_index_map_single_component, _mesh);
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
        mesh.getElements(), dof_table, _local_assemblers,
        mesh.isAxiallySymmetric(), integration_order, _process_data);

    _secondary_variables.addSecondaryVariable(
        "sigma",
        makeExtrapolator(
            MathLib::KelvinVector::KelvinVectorType<
                DisplacementDim>::RowsAtCompileTime,
            getExtrapolator(), _local_assemblers,
            &ThermoMechanicsLocalAssemblerInterface::getIntPtSigma));

    _secondary_variables.addSecondaryVariable(
        "epsilon",
        makeExtrapolator(
            MathLib::KelvinVector::KelvinVectorType<
                DisplacementDim>::RowsAtCompileTime,
            getExtrapolator(), _local_assemblers,
            &ThermoMechanicsLocalAssemblerInterface::getIntPtEpsilon));

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
        if (ip_meta_data.n_components != mesh_property.getNumberOfComponents())
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

    // Initialize local assemblers after all variables have been set.
    GlobalExecutor::executeMemberOnDereferenced(
        &LocalAssemblerInterface::initialize, _local_assemblers,
        *_local_to_global_index_map);
}

template <int DisplacementDim>
void ThermoMechanicsProcess<DisplacementDim>::initializeBoundaryConditions()
{
    if (_use_monolithic_scheme)
    {
        const int process_id_of_thermomechanics = 0;
        initializeProcessBoundaryConditionsAndSourceTerms(
            *_local_to_global_index_map, process_id_of_thermomechanics);
        return;
    }

    // Staggered scheme:
    // for the equations of heat conduction
    initializeProcessBoundaryConditionsAndSourceTerms(
        *_local_to_global_index_map_single_component,
        _process_data.heat_conduction_process_id);

    // for the equations of deformation.
    initializeProcessBoundaryConditionsAndSourceTerms(
        *_local_to_global_index_map, _process_data.mechanics_process_id);
}

template <int DisplacementDim>
void ThermoMechanicsProcess<DisplacementDim>::assembleConcreteProcess(
    const double t, GlobalVector const& x, GlobalMatrix& M, GlobalMatrix& K,
    GlobalVector& b)
{
    DBUG("Assemble ThermoMechanicsProcess.");

    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_table = {std::ref(*_local_to_global_index_map)};
    const int process_id =
        _use_monolithic_scheme ? 0 : _coupled_solutions->process_id;
    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    // Call global assembler for each local assembly item.
    GlobalExecutor::executeSelectedMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assemble, _local_assemblers,
        pv.getActiveElementIDs(), dof_table, t, x, M, K, b, _coupled_solutions);
}

template <int DisplacementDim>
void ThermoMechanicsProcess<DisplacementDim>::
    assembleWithJacobianConcreteProcess(const double t, GlobalVector const& x,
                                        GlobalVector const& xdot,
                                        const double dxdot_dx,
                                        const double dx_dx, GlobalMatrix& M,
                                        GlobalMatrix& K, GlobalVector& b,
                                        GlobalMatrix& Jac)
{
    DBUG("AssembleJacobian ThermoMechanicsProcess.");

    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_tables;
    // For the monolithic scheme
    if (_use_monolithic_scheme)
    {
        DBUG(
            "Assemble the Jacobian of ThermoMechanics for the monolithic"
            " scheme.");
        dof_tables.emplace_back(*_local_to_global_index_map);
    }
    else
    {
        // For the staggered scheme
        if (_coupled_solutions->process_id ==
            _process_data.heat_conduction_process_id)
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
        if (_process_data.heat_conduction_process_id ==
            0)  // First: the heat conduction process
        {
            dof_tables.emplace_back(
                *_local_to_global_index_map_single_component);
            dof_tables.emplace_back(*_local_to_global_index_map);
        }
        else  // vice versa
        {
            dof_tables.emplace_back(*_local_to_global_index_map);
            dof_tables.emplace_back(
                *_local_to_global_index_map_single_component);
        }

        setCoupledSolutionsOfPreviousTimeStep();
    }

    const int process_id =
        _use_monolithic_scheme ? 0 : _coupled_solutions->process_id;
    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    GlobalExecutor::executeSelectedMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assembleWithJacobian,
        _local_assemblers, pv.getActiveElementIDs(), dof_tables, t, x, xdot,
        dxdot_dx, dx_dx, M, K, b, Jac, _coupled_solutions);

    // TODO (naumov): Refactor the copy rhs part. This is copy from HM.
    auto copyRhs = [&](int const variable_id, auto& output_vector) {
        if (_use_monolithic_scheme)
        {
            transformVariableFromGlobalVector(b, variable_id, dof_tables[0],
                                              output_vector,
                                              std::negate<double>());
        }
        else
        {
            transformVariableFromGlobalVector(
                b, 0, dof_tables[_coupled_solutions->process_id], output_vector,
                std::negate<double>());
        }
    };
    if (_use_monolithic_scheme ||
        _coupled_solutions->process_id ==
            _process_data.heat_conduction_process_id)
    {
        copyRhs(0, *_heat_flux);
    }
    if (_use_monolithic_scheme ||
        _coupled_solutions->process_id == _process_data.mechanics_process_id)
    {
        copyRhs(1, *_nodal_forces);
    }
}

template <int DisplacementDim>
void ThermoMechanicsProcess<DisplacementDim>::preTimestepConcreteProcess(
    GlobalVector const& x, double const t, double const dt,
    const int process_id)
{
    DBUG("PreTimestep ThermoMechanicsProcess.");

    _process_data.dt = dt;
    _process_data.t = t;

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    assert(process_id < 2);

    if (process_id == _process_data.mechanics_process_id)
    {
        GlobalExecutor::executeSelectedMemberOnDereferenced(
            &ThermoMechanicsLocalAssemblerInterface::preTimestep,
            _local_assemblers, pv.getActiveElementIDs(),
            *_local_to_global_index_map, x, t, dt);
        return;
    }

    // For the staggered scheme.
    if (!_previous_T)
    {
        _previous_T = MathLib::MatrixVectorTraits<GlobalVector>::newInstance(x);
    }
    else
    {
        auto& x0 = *_previous_T;
        MathLib::LinAlg::copy(x, x0);
    }

    auto& x0 = *_previous_T;
    MathLib::LinAlg::setLocalAccessibleVector(x0);
}

template <int DisplacementDim>
void ThermoMechanicsProcess<DisplacementDim>::postTimestepConcreteProcess(
    GlobalVector const& x, const double /*t*/, const double /*delta_t*/,
    int const process_id)
{
    if (process_id != _process_data.mechanics_process_id)
    {
        return;
    }

    DBUG("PostTimestep ThermoMechanicsProcess.");

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    GlobalExecutor::executeSelectedMemberOnDereferenced(
        &ThermoMechanicsLocalAssemblerInterface::postTimestep,
        _local_assemblers, pv.getActiveElementIDs(),
        *_local_to_global_index_map, x);
}

template <int DisplacementDim>
void ThermoMechanicsProcess<
    DisplacementDim>::setCoupledSolutionsOfPreviousTimeStep()
{
    _coupled_solutions->coupled_xs_t0.resize(1);
    _coupled_solutions->coupled_xs_t0[0] = _previous_T.get();
}

template <int DisplacementDim>
NumLib::LocalToGlobalIndexMap const&
ThermoMechanicsProcess<DisplacementDim>::getDOFTable(const int process_id) const
{
    if (_process_data.mechanics_process_id == process_id)
    {
        return *_local_to_global_index_map;
    }

    // For the equation of pressure
    return *_local_to_global_index_map_single_component;
}

template class ThermoMechanicsProcess<2>;
template class ThermoMechanicsProcess<3>;

}  // namespace ThermoMechanics
}  // namespace ProcessLib
