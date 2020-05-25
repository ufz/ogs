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
      process_data_(std::move(process_data)),
      surfaceflux_(std::move(surfaceflux))
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
    if (use_monolithic_scheme_)
    {
        for (auto pv_iter = std::next(process_variables_[process_id].begin());
             pv_iter != process_variables_[process_id].end();
             ++pv_iter)
        {
            transport_process_variables.push_back(*pv_iter);
        }
    }
    else
    {
        for (auto pv_iter = std::next(process_variables_.begin());
             pv_iter != process_variables_.end();
             ++pv_iter)
        {
            transport_process_variables.push_back((*pv_iter)[0]);
        }

        xs_previous_timestep_.resize(process_variables_.size());
    }

    ProcessLib::createLocalAssemblers<LocalAssemblerData>(
        mesh.getDimension(), mesh.getElements(), dof_table,
        pv.getShapeFunctionOrder(), local_assemblers_,
        mesh.isAxiallySymmetric(), integration_order, process_data_,
        transport_process_variables);

    secondary_variables_.addSecondaryVariable(
        "darcy_velocity",
        makeExtrapolator(
            mesh.getDimension(), getExtrapolator(), local_assemblers_,
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
    if (use_monolithic_scheme_)
    {
        dof_tables.push_back(std::ref(*local_to_global_index_map_));
    }
    else
    {
        setCoupledSolutionsOfPreviousTimeStep();
        std::generate_n(
            std::back_inserter(dof_tables), process_variables_.size(),
            [&]() { return std::ref(*local_to_global_index_map_); });
    }
    // Call global assembler for each local assembly item.
    GlobalExecutor::executeSelectedMemberDereferenced(
        global_assembler_, &VectorMatrixAssembler::assemble, local_assemblers_,
        pv.getActiveElementIDs(), dof_tables, t, dt, x, xdot, process_id, M, K,
        b, coupled_solutions_);
}

void ComponentTransportProcess::setCoupledSolutionsOfPreviousTimeStep()
{
    unsigned const number_of_coupled_solutions =
        coupled_solutions_->coupled_xs.size();
    coupled_solutions_->coupled_xs_t0.clear();
    coupled_solutions_->coupled_xs_t0.reserve(number_of_coupled_solutions);
    for (unsigned i = 0; i < number_of_coupled_solutions; ++i)
    {
        auto const& x_t0 = xs_previous_timestep_[i];
        coupled_solutions_->coupled_xs_t0.emplace_back(x_t0.get());
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
        dof_table = {std::ref(*local_to_global_index_map_)};
    // Call global assembler for each local assembly item.
    GlobalExecutor::executeSelectedMemberDereferenced(
        global_assembler_, &VectorMatrixAssembler::assembleWithJacobian,
        local_assemblers_, pv.getActiveElementIDs(), dof_table, t, dt, x, xdot,
        dxdot_dx, dx_dx, process_id, M, K, b, Jac, coupled_solutions_);
}

Eigen::Vector3d ComponentTransportProcess::getFlux(
    std::size_t const element_id,
    MathLib::Point3d const& p,
    double const t,
    std::vector<GlobalVector*> const& x) const
{
    std::vector<GlobalIndexType> indices_cache;
    auto const r_c_indices = NumLib::getRowColumnIndices(
        element_id, *local_to_global_index_map_, indices_cache);

    std::vector<std::vector<GlobalIndexType>> indices_of_all_coupled_processes{
        x.size(), r_c_indices.rows};
    auto const local_xs =
        getCoupledLocalSolutions(x, indices_of_all_coupled_processes);

    return local_assemblers_[element_id]->getFlux(p, t, local_xs);
}

void ComponentTransportProcess::
    setCoupledTermForTheStaggeredSchemeToLocalAssemblers(int const process_id)
{
    DBUG("Set the coupled term for the staggered scheme to local assembers.");

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];
    GlobalExecutor::executeSelectedMemberOnDereferenced(
        &ComponentTransportLocalAssemblerInterface::
            setStaggeredCoupledSolutions,
        local_assemblers_, pv.getActiveElementIDs(), coupled_solutions_);
}

void ComponentTransportProcess::preTimestepConcreteProcess(
    std::vector<GlobalVector*> const& x, const double /*t*/,
    const double /*delta_t*/, int const process_id)
{
    if (use_monolithic_scheme_)
    {
        return;
    }

    if (!xs_previous_timestep_[process_id])
    {
        xs_previous_timestep_[process_id] =
            MathLib::MatrixVectorTraits<GlobalVector>::newInstance(
                *x[process_id]);
    }
    else
    {
        auto& x0 = *xs_previous_timestep_[process_id];
        MathLib::LinAlg::copy(*x[process_id], x0);
    }

    auto& x0 = *xs_previous_timestep_[process_id];
    MathLib::LinAlg::setLocalAccessibleVector(x0);
}

void ComponentTransportProcess::postTimestepConcreteProcess(
    std::vector<GlobalVector*> const& x,
    const double t,
    const double /*delta_t*/,
    int const process_id)
{
    // For the monolithic scheme, process_id is always zero.
    if (use_monolithic_scheme_ && process_id != 0)
    {
        OGS_FATAL(
            "The condition of process_id = 0 must be satisfied for "
            "monolithic ComponentTransportProcess, which is a single process.");
    }
    if (!use_monolithic_scheme_ && process_id != 0)
    {
        DBUG(
            "This is the transport part of the staggered "
            "ComponentTransportProcess.");
        return;
    }
    if (!surfaceflux_)  // computing the surfaceflux is optional
    {
        return;
    }

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    surfaceflux_->integrate(x, t, *this, process_id, integration_order_, mesh_,
                            pv.getActiveElementIDs());
    surfaceflux_->save(t);
}

}  // namespace ComponentTransport
}  // namespace ProcessLib
