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
    std::unique_ptr<ProcessLib::SurfaceFluxData>&& surfaceflux)
    : Process(std::move(name), mesh, std::move(jacobian_assembler), parameters,
              integration_order, std::move(process_variables),
              std::move(secondary_variables), use_monolithic_scheme),
      _process_data(std::move(process_data)),
      _surfaceflux(std::move(surfaceflux))
{
}

void ComponentTransportProcess::initializeConcreteProcess(
    NumLib::LocalToGlobalIndexMap const& dof_table,
    MeshLib::Mesh const& mesh,
    unsigned const integration_order)
{
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

        _xs_previous_timestep.resize(_process_variables.size());
    }

    ProcessLib::createLocalAssemblers<LocalAssemblerData>(
        mesh.getDimension(), mesh.getElements(), dof_table,
        pv.getShapeFunctionOrder(), _local_assemblers,
        mesh.isAxiallySymmetric(), integration_order, _process_data,
        transport_process_variables);

    _secondary_variables.addSecondaryVariable(
        "darcy_velocity",
        makeExtrapolator(
            mesh.getDimension(), getExtrapolator(), _local_assemblers,
            &ComponentTransportLocalAssemblerInterface::getIntPtDarcyVelocity));
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
        setCoupledSolutionsOfPreviousTimeStep();
        std::generate_n(
            std::back_inserter(dof_tables), _process_variables.size(),
            [&]() { return std::ref(*_local_to_global_index_map); });
    }
    // Call global assembler for each local assembly item.
    GlobalExecutor::executeSelectedMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assemble, _local_assemblers,
        pv.getActiveElementIDs(), dof_tables, t, dt, x, xdot, process_id, M, K,
        b, _coupled_solutions);
}

void ComponentTransportProcess::setCoupledSolutionsOfPreviousTimeStep()
{
    unsigned const number_of_coupled_solutions =
        _coupled_solutions->coupled_xs.size();
    _coupled_solutions->coupled_xs_t0.clear();
    _coupled_solutions->coupled_xs_t0.reserve(number_of_coupled_solutions);
    for (unsigned i = 0; i < number_of_coupled_solutions; ++i)
    {
        auto const& x_t0 = _xs_previous_timestep[i];
        _coupled_solutions->coupled_xs_t0.emplace_back(x_t0.get());
    }
}

void ComponentTransportProcess::assembleWithJacobianConcreteProcess(
    const double t, double const dt, std::vector<GlobalVector*> const& x,
    GlobalVector const& xdot, const double dxdot_dx, const double dx_dx,
    int const process_id, GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b,
    GlobalMatrix& Jac)
{
    DBUG("AssembleWithJacobian ComponentTransportProcess.");

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];
    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_table = {std::ref(*_local_to_global_index_map)};
    // Call global assembler for each local assembly item.
    GlobalExecutor::executeSelectedMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assembleWithJacobian,
        _local_assemblers, pv.getActiveElementIDs(), dof_table, t, dt, x, xdot,
        dxdot_dx, dx_dx, process_id, M, K, b, Jac, _coupled_solutions);
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

std::vector<GlobalVector>
ComponentTransportProcess::interpolateProcessSolutions(
    std::vector<GlobalVector*> const& x)
{
    // Result is for each process a vector of integration point values for each
    // element stored consecutively.
    auto interpolateNodalValuesToIntegrationPoints =
        [this](std::size_t mesh_item_id,
               LocalAssemblerInterface& local_assembler,
               std::vector<
                   std::reference_wrapper<NumLib::LocalToGlobalIndexMap>> const&
                   dof_tables,
               std::vector<GlobalVector*> const& x,
               std::vector<std::vector<double>>& int_pt_process_solutions) {
            for (unsigned process_id = 0; process_id < x.size(); ++process_id)
            {
                auto const& dof_table = dof_tables[process_id].get();
                auto const indices =
                    NumLib::getIndices(mesh_item_id, dof_table);
                auto const local_x = x[process_id]->get(indices);

                // interpolated_values in size of integration points
                std::vector<double> const interpolated_values =
                    local_assembler.interpolateNodalValuesToIntegrationPoints(
                        local_x);

                // For each element (mesh_item_id) concatenate the integration
                // point values.
                auto& int_pt_process_solution =
                    int_pt_process_solutions[process_id];
                int_pt_process_solution.insert(int_pt_process_solution.end(),
                                               interpolated_values.begin(),
                                               interpolated_values.end());
            }
        };

    // Same dof table for each primary variable.
    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_tables;
    dof_tables.reserve(x.size());
    std::generate_n(std::back_inserter(dof_tables), x.size(),
                    [&]() { return std::ref(*_local_to_global_index_map); });

    std::vector<std::vector<double>> int_pt_process_solutions(x.size());

    GlobalExecutor::executeDereferenced(
        interpolateNodalValuesToIntegrationPoints, _local_assemblers,
        dof_tables, x, int_pt_process_solutions);

    // For each process copy the (elements x integration_points_per_element)
    // vector to the result `int_pt_x`.
    std::vector<GlobalVector> int_pt_x;
    int_pt_x.reserve(x.size());
    for (auto const& int_pt_process_solution : int_pt_process_solutions)
    {
        GlobalIndexType const size = int_pt_process_solution.size();
        // New vector of size (elements x integration_points_per_element)
        GlobalVector int_pt_process_solution_vec(size);

        // Copy one by one.
        for (GlobalIndexType i = 0; i < size; ++i)
        {
            int_pt_process_solution_vec.set(i, int_pt_process_solution[i]);
        }
        int_pt_x.push_back(int_pt_process_solution_vec);
    }

    return int_pt_x;
}

void ComponentTransportProcess::extrapolateIntegrationPointValuesToNodes(
    const double t, std::vector<GlobalVector*> const& int_pt_x,
    std::vector<GlobalVector*>& x)
{
    auto& extrapolator = getExtrapolator();
    auto const extrapolatables =
        NumLib::makeExtrapolatable(_local_assemblers,
                                   &ComponentTransportLocalAssemblerInterface::
                                       getInterpolatedLocalSolution);

    for (unsigned transport_process_id = 0;
         transport_process_id < int_pt_x.size();
         ++transport_process_id)
    {
        auto const& pv = _process_variables[transport_process_id][0].get();
        auto const& int_pt_C = int_pt_x[transport_process_id];

        extrapolator.extrapolate(pv.getNumberOfComponents(), extrapolatables, t,
                                 {int_pt_C},
                                 {_local_to_global_index_map.get()});

        auto const& nodal_values = extrapolator.getNodalValues();
        *x[transport_process_id + 1] = nodal_values;
    }
}

void ComponentTransportProcess::preTimestepConcreteProcess(
    std::vector<GlobalVector*> const& x, const double /*t*/,
    const double /*delta_t*/, int const process_id)
{
    if (_use_monolithic_scheme)
    {
        return;
    }

    if (!_xs_previous_timestep[process_id])
    {
        _xs_previous_timestep[process_id] =
            MathLib::MatrixVectorTraits<GlobalVector>::newInstance(
                *x[process_id]);
    }
    else
    {
        auto& x0 = *_xs_previous_timestep[process_id];
        MathLib::LinAlg::copy(*x[process_id], x0);
    }

    auto& x0 = *_xs_previous_timestep[process_id];
    MathLib::LinAlg::setLocalAccessibleVector(x0);
}

void ComponentTransportProcess::postTimestepConcreteProcess(
    std::vector<GlobalVector*> const& x,
    const double t,
    const double /*delta_t*/,
    int const process_id)
{
    // For the monolithic scheme, process_id is always zero.
    if (_use_monolithic_scheme && process_id != 0)
    {
        OGS_FATAL(
            "The condition of process_id = 0 must be satisfied for "
            "monolithic ComponentTransportProcess, which is a single process.");
    }
    if (!_use_monolithic_scheme && process_id != 0)
    {
        DBUG(
            "This is the transport part of the staggered "
            "ComponentTransportProcess.");
        return;
    }
    if (!_surfaceflux)  // computing the surfaceflux is optional
    {
        return;
    }

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    _surfaceflux->integrate(x, t, *this, process_id, _integration_order, _mesh,
                            pv.getActiveElementIDs());
    _surfaceflux->save(t);
}

}  // namespace ComponentTransport
}  // namespace ProcessLib
