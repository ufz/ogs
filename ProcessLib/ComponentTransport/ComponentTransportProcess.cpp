/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ComponentTransportProcess.h"

#include <cassert>

#include "BaseLib/RunTime.h"
#include "ChemistryLib/ChemicalSolverInterface.h"
#include "ProcessLib/SurfaceFlux/SurfaceFlux.h"
#include "ProcessLib/SurfaceFlux/SurfaceFluxData.h"
#include "ProcessLib/Utils/CreateLocalAssemblers.h"

namespace ProcessLib
{
namespace ComponentTransport
{
ComponentTransportProcess::ComponentTransportProcess(
    std::string name,
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    unsigned const integration_order,
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
        process_variables,
    ComponentTransportProcessData&& process_data,
    SecondaryVariableCollection&& secondary_variables,
    bool const use_monolithic_scheme,
    std::unique_ptr<ProcessLib::SurfaceFluxData>&& surfaceflux,
    std::unique_ptr<ChemistryLib::ChemicalSolverInterface>&&
        chemical_solver_interface)
    : Process(std::move(name), mesh, std::move(jacobian_assembler), parameters,
              integration_order, std::move(process_variables),
              std::move(secondary_variables), use_monolithic_scheme),
      _process_data(std::move(process_data)),
      _surfaceflux(std::move(surfaceflux)),
      _chemical_solver_interface(std::move(chemical_solver_interface))
{
}

void ComponentTransportProcess::initializeConcreteProcess(
    NumLib::LocalToGlobalIndexMap const& dof_table,
    MeshLib::Mesh const& mesh,
    unsigned const integration_order)
{
    _process_data.mesh_prop_velocity = MeshLib::getOrCreateMeshProperty<double>(
        const_cast<MeshLib::Mesh&>(mesh), "velocity",
        MeshLib::MeshItemType::Cell, mesh.getDimension());

    const int process_id = 0;
    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    std::vector<std::reference_wrapper<ProcessLib::ProcessVariable>>
        transport_process_variables;
    if (_use_monolithic_scheme)
    {
        for (auto pv_iter = std::next(_process_variables[process_id].begin());
             pv_iter != _process_variables[process_id].end();
             ++pv_iter)
        {
            transport_process_variables.push_back(*pv_iter);
        }
    }
    else
    {
        for (auto pv_iter = std::next(_process_variables.begin());
             pv_iter != _process_variables.end();
             ++pv_iter)
        {
            transport_process_variables.push_back((*pv_iter)[0]);
        }
    }

    ProcessLib::createLocalAssemblers<LocalAssemblerData>(
        mesh.getDimension(), mesh.getElements(), dof_table,
        pv.getShapeFunctionOrder(), _local_assemblers,
        mesh.isAxiallySymmetric(), integration_order, _process_data,
        transport_process_variables);

    if (_chemical_solver_interface)
    {
        _chemical_solver_interface->initialize();
    }

    _secondary_variables.addSecondaryVariable(
        "darcy_velocity",
        makeExtrapolator(
            mesh.getDimension(), getExtrapolator(), _local_assemblers,
            &ComponentTransportLocalAssemblerInterface::getIntPtDarcyVelocity));
}

void ComponentTransportProcess::setInitialConditionsConcreteProcess(
    std::vector<GlobalVector*>& x, double const t, int const process_id)
{
    if (!_chemical_solver_interface)
    {
        return;
    }

    if (process_id != static_cast<int>(x.size() - 1))
    {
        return;
    }

    BaseLib::RunTime time_phreeqc;
    time_phreeqc.start();

    _chemical_solver_interface->executeInitialCalculation(
        interpolateNodalValuesToIntegrationPoints(x));

    extrapolateIntegrationPointValuesToNodes(
        t, _chemical_solver_interface->getIntPtProcessSolutions(), x);

    INFO("[time] Phreeqc took {:g} s.", time_phreeqc.elapsed());
}

void ComponentTransportProcess::assembleConcreteProcess(
    const double t, double const dt, std::vector<GlobalVector*> const& x,
    std::vector<GlobalVector*> const& xdot, int const process_id,
    GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b)
{
    DBUG("Assemble ComponentTransportProcess.");

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_tables;
    if (_use_monolithic_scheme)
    {
        dof_tables.push_back(std::ref(*_local_to_global_index_map));
    }
    else
    {
        std::generate_n(
            std::back_inserter(dof_tables), _process_variables.size(),
            [&]() { return std::ref(*_local_to_global_index_map); });
    }
    // Call global assembler for each local assembly item.
    GlobalExecutor::executeSelectedMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assemble, _local_assemblers,
        pv.getActiveElementIDs(), dof_tables, t, dt, x, xdot, process_id, M, K,
        b);
}

void ComponentTransportProcess::assembleWithJacobianConcreteProcess(
    const double t, double const dt, std::vector<GlobalVector*> const& x,
    std::vector<GlobalVector*> const& xdot, const double dxdot_dx,
    const double dx_dx, int const process_id, GlobalMatrix& M, GlobalMatrix& K,
    GlobalVector& b, GlobalMatrix& Jac)
{
    DBUG("AssembleWithJacobian ComponentTransportProcess.");

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];
    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_table = {std::ref(*_local_to_global_index_map)};
    // Call global assembler for each local assembly item.
    GlobalExecutor::executeSelectedMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assembleWithJacobian,
        _local_assemblers, pv.getActiveElementIDs(), dof_table, t, dt, x, xdot,
        dxdot_dx, dx_dx, process_id, M, K, b, Jac);
}

Eigen::Vector3d ComponentTransportProcess::getFlux(
    std::size_t const element_id,
    MathLib::Point3d const& p,
    double const t,
    std::vector<GlobalVector*> const& x) const
{
    std::vector<GlobalIndexType> indices_cache;
    auto const r_c_indices = NumLib::getRowColumnIndices(
        element_id, *_local_to_global_index_map, indices_cache);

    std::vector<std::vector<GlobalIndexType>> indices_of_all_coupled_processes{
        x.size(), r_c_indices.rows};
    auto const local_xs =
        getCoupledLocalSolutions(x, indices_of_all_coupled_processes);

    return _local_assemblers[element_id]->getFlux(p, t, local_xs);
}

void ComponentTransportProcess::
    setCoupledTermForTheStaggeredSchemeToLocalAssemblers(int const process_id)
{
    DBUG("Set the coupled term for the staggered scheme to local assembers.");

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];
    GlobalExecutor::executeSelectedMemberOnDereferenced(
        &ComponentTransportLocalAssemblerInterface::
            setStaggeredCoupledSolutions,
        _local_assemblers, pv.getActiveElementIDs(), _coupled_solutions);
}

void ComponentTransportProcess::solveReactionEquation(
    std::vector<GlobalVector*>& x, double const t, double const dt)
{
    if (!_chemical_solver_interface)
    {
        return;
    }

    // Sequential non-iterative approach applied here to perform water
    // chemistry calculation followed by resolving component transport
    // process.
    // TODO: move into a global loop to consider both mass balance over
    // space and localized chemical equilibrium between solutes.
    BaseLib::RunTime time_phreeqc;
    time_phreeqc.start();

    _chemical_solver_interface->doWaterChemistryCalculation(
        interpolateNodalValuesToIntegrationPoints(x), dt);

    extrapolateIntegrationPointValuesToNodes(
        t, _chemical_solver_interface->getIntPtProcessSolutions(), x);

    INFO("[time] Phreeqc took {:g} s.", time_phreeqc.elapsed());
}

std::vector<GlobalVector>
ComponentTransportProcess::interpolateNodalValuesToIntegrationPoints(
    std::vector<GlobalVector*> const& nodal_values_vectors) const
{
    // Result is for each process a vector of integration point values for each
    // element stored consecutively.
    auto interpolateNodalValuesToIntegrationPoints =
        [](std::size_t mesh_item_id,
           LocalAssemblerInterface& local_assembler,
           std::vector<
               std::reference_wrapper<NumLib::LocalToGlobalIndexMap>> const&
               dof_tables,
           std::vector<GlobalVector*> const& x,
           std::vector<std::vector<double>>& int_pt_x) {
            for (unsigned process_id = 0; process_id < x.size(); ++process_id)
            {
                auto const& dof_table = dof_tables[process_id].get();
                auto const indices =
                    NumLib::getIndices(mesh_item_id, dof_table);
                auto const local_x = x[process_id]->get(indices);

                std::vector<double> const interpolated_values =
                    local_assembler.interpolateNodalValuesToIntegrationPoints(
                        local_x);

                // For each element (mesh_item_id) concatenate the integration
                // point values.
                int_pt_x[process_id].insert(int_pt_x[process_id].end(),
                                            interpolated_values.begin(),
                                            interpolated_values.end());
            }
        };

    // Same dof table for each primary variable.
    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_tables;
    dof_tables.reserve(nodal_values_vectors.size());
    std::generate_n(std::back_inserter(dof_tables), nodal_values_vectors.size(),
                    [&]() { return std::ref(*_local_to_global_index_map); });

    std::vector<std::vector<double>> integration_point_values_vecs(
        nodal_values_vectors.size());

    GlobalExecutor::executeDereferenced(
        interpolateNodalValuesToIntegrationPoints, _local_assemblers,
        dof_tables, nodal_values_vectors, integration_point_values_vecs);

    // Convert from std::vector<std::vector<double>> to
    // std::vector<GlobalVector>
    std::vector<GlobalVector> integration_point_values_vectors;
    integration_point_values_vectors.reserve(nodal_values_vectors.size());
    for (auto const& integration_point_values_vec :
         integration_point_values_vecs)
    {
        GlobalIndexType const size = integration_point_values_vec.size();
        // New vector of size (elements x integration_points_per_element)
        GlobalVector integration_point_values_vector(size);

        // Copy one by one.
        for (GlobalIndexType i = 0; i < size; ++i)
        {
            integration_point_values_vector.set(
                i, integration_point_values_vec[i]);
        }
        integration_point_values_vectors.push_back(
            std::move(integration_point_values_vector));
    }

    return integration_point_values_vectors;
}

void ComponentTransportProcess::extrapolateIntegrationPointValuesToNodes(
    const double t,
    std::vector<GlobalVector*> const& integration_point_values_vectors,
    std::vector<GlobalVector*>& nodal_values_vectors)
{
    auto& extrapolator = getExtrapolator();
    auto const extrapolatables =
        NumLib::makeExtrapolatable(_local_assemblers,
                                   &ComponentTransportLocalAssemblerInterface::
                                       getInterpolatedLocalSolution);

    for (unsigned transport_process_id = 0;
         transport_process_id < integration_point_values_vectors.size();
         ++transport_process_id)
    {
        auto const& pv = _process_variables[transport_process_id + 1][0].get();
        auto const& int_pt_C =
            integration_point_values_vectors[transport_process_id];

        extrapolator.extrapolate(pv.getNumberOfGlobalComponents(),
                                 extrapolatables, t, {int_pt_C},
                                 {_local_to_global_index_map.get()});

        auto const& nodal_values = extrapolator.getNodalValues();
        MathLib::LinAlg::copy(nodal_values,
                              *nodal_values_vectors[transport_process_id + 1]);
    }
}

void ComponentTransportProcess::computeSecondaryVariableConcrete(
    double const t,
    double const dt,
    std::vector<GlobalVector*> const& x,
    GlobalVector const& x_dot,
    int const process_id)
{
    if (process_id != 0)
    {
        return;
    }

    std::vector<NumLib::LocalToGlobalIndexMap const*> dof_tables;
    dof_tables.reserve(x.size());
    std::generate_n(std::back_inserter(dof_tables), x.size(),
                    [&]() { return _local_to_global_index_map.get(); });

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];
    GlobalExecutor::executeSelectedMemberOnDereferenced(
        &ComponentTransportLocalAssemblerInterface::computeSecondaryVariable,
        _local_assemblers, pv.getActiveElementIDs(), dof_tables, t, dt, x,
        x_dot, process_id);
}

void ComponentTransportProcess::postTimestepConcreteProcess(
    std::vector<GlobalVector*> const& x,
    const double t,
    const double dt,
    int const process_id)
{
    if (process_id != 0)
    {
        return;
    }

    std::vector<NumLib::LocalToGlobalIndexMap const*> dof_tables;
    dof_tables.reserve(x.size());
    std::generate_n(std::back_inserter(dof_tables), x.size(),
                    [&]() { return _local_to_global_index_map.get(); });

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];
    GlobalExecutor::executeSelectedMemberOnDereferenced(
        &ComponentTransportLocalAssemblerInterface::postTimestep,
        _local_assemblers, pv.getActiveElementIDs(), dof_tables, x, t, dt);

    if (!_surfaceflux)  // computing the surfaceflux is optional
    {
        return;
    }
    _surfaceflux->integrate(x, t, *this, process_id, _integration_order, _mesh,
                            pv.getActiveElementIDs());
}

}  // namespace ComponentTransport
}  // namespace ProcessLib
