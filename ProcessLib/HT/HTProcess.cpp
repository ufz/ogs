/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "HTProcess.h"

#include <cassert>

#include "NumLib/DOF/DOFTableUtil.h"
#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "ProcessLib/SurfaceFlux/SurfaceFluxData.h"
#include "ProcessLib/Utils/CreateLocalAssemblers.h"

#include "MonolithicHTFEM.h"
#include "StaggeredHTFEM.h"

namespace ProcessLib
{
namespace HT
{
HTProcess::HTProcess(
    std::string name,
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    unsigned const integration_order,
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
        process_variables,
    HTProcessData&& process_data,
    SecondaryVariableCollection&& secondary_variables,
    bool const use_monolithic_scheme,
    std::unique_ptr<ProcessLib::SurfaceFluxData>&& surfaceflux,
    const int heat_transport_process_id,
    const int hydraulic_process_id)
    : Process(std::move(name), mesh, std::move(jacobian_assembler), parameters,
              integration_order, std::move(process_variables),
              std::move(secondary_variables), use_monolithic_scheme),
      process_data_(std::move(process_data)),
      surfaceflux_(std::move(surfaceflux)),
      heat_transport_process_id_(heat_transport_process_id),
      hydraulic_process_id_(hydraulic_process_id)
{
}

void HTProcess::initializeConcreteProcess(
    NumLib::LocalToGlobalIndexMap const& dof_table,
    MeshLib::Mesh const& mesh,
    unsigned const integration_order)
{
    // For the staggered scheme, both processes are assumed to use the same
    // element order. Therefore the order of shape function can be fetched from
    // any set of the sets of process variables of the coupled processes. Here,
    // we take the one from the first process by setting process_id = 0.
    const int process_id = 0;
    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    if (use_monolithic_scheme_)
    {
        ProcessLib::createLocalAssemblers<MonolithicHTFEM>(
            mesh.getDimension(), mesh.getElements(), dof_table,
            pv.getShapeFunctionOrder(), local_assemblers_,
            mesh.isAxiallySymmetric(), integration_order, process_data_);
    }
    else
    {
        ProcessLib::createLocalAssemblers<StaggeredHTFEM>(
            mesh.getDimension(), mesh.getElements(), dof_table,
            pv.getShapeFunctionOrder(), local_assemblers_,
            mesh.isAxiallySymmetric(), integration_order, process_data_,
            heat_transport_process_id_, hydraulic_process_id_);
    }

    secondary_variables_.addSecondaryVariable(
        "darcy_velocity",
        makeExtrapolator(mesh.getDimension(), getExtrapolator(),
                         local_assemblers_,
                         &HTLocalAssemblerInterface::getIntPtDarcyVelocity));
}

void HTProcess::assembleConcreteProcess(const double t, double const dt,
                                        std::vector<GlobalVector*> const& x,
                                        std::vector<GlobalVector*> const& xdot,
                                        int const process_id, GlobalMatrix& M,
                                        GlobalMatrix& K, GlobalVector& b)
{
    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_tables;
    if (use_monolithic_scheme_)
    {
        DBUG("Assemble HTProcess.");
        dof_tables.emplace_back(*local_to_global_index_map_);
    }
    else
    {
        if (process_id == heat_transport_process_id_)
        {
            DBUG(
                "Assemble the equations of heat transport process within "
                "HTProcess.");
        }
        else
        {
            DBUG(
                "Assemble the equations of single phase fully saturated "
                "fluid flow process within HTProcess.");
        }
        setCoupledSolutionsOfPreviousTimeStep();
        dof_tables.emplace_back(*local_to_global_index_map_);
        dof_tables.emplace_back(*local_to_global_index_map_);
    }

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];
    // Call global assembler for each local assembly item.
    GlobalExecutor::executeSelectedMemberDereferenced(
        global_assembler_, &VectorMatrixAssembler::assemble, local_assemblers_,
        pv.getActiveElementIDs(), dof_tables, t, dt, x, xdot, process_id, M, K,
        b, coupled_solutions_);
}

void HTProcess::assembleWithJacobianConcreteProcess(
    const double t, double const dt, std::vector<GlobalVector*> const& x,
    GlobalVector const& xdot, const double dxdot_dx, const double dx_dx,
    int const process_id, GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b,
    GlobalMatrix& Jac)
{
    DBUG("AssembleWithJacobian HTProcess.");

    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_tables;
    if (!use_monolithic_scheme_)
    {
        setCoupledSolutionsOfPreviousTimeStep();
        dof_tables.emplace_back(std::ref(*local_to_global_index_map_));
    }
    else
    {
        dof_tables.emplace_back(std::ref(*local_to_global_index_map_));
        dof_tables.emplace_back(std::ref(*local_to_global_index_map_));
    }

    // Call global assembler for each local assembly item.
    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];
    GlobalExecutor::executeSelectedMemberDereferenced(
        global_assembler_, &VectorMatrixAssembler::assembleWithJacobian,
        local_assemblers_, pv.getActiveElementIDs(), dof_tables, t, dt, x, xdot,
        dxdot_dx, dx_dx, process_id, M, K, b, Jac, coupled_solutions_);
}

void HTProcess::preTimestepConcreteProcess(std::vector<GlobalVector*> const& x,
                                           const double /*t*/,
                                           const double /*delta_t*/,
                                           const int process_id)
{
    assert(process_id < 2);

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

void HTProcess::setCoupledTermForTheStaggeredSchemeToLocalAssemblers(
    int const process_id)
{
    DBUG("Set the coupled term for the staggered scheme to local assembers.");

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];
    GlobalExecutor::executeSelectedMemberOnDereferenced(
        &HTLocalAssemblerInterface::setStaggeredCoupledSolutions,
        local_assemblers_, pv.getActiveElementIDs(), coupled_solutions_);
}

std::tuple<NumLib::LocalToGlobalIndexMap*, bool>
HTProcess::getDOFTableForExtrapolatorData() const
{
    if (!use_monolithic_scheme_)
    {
        // For single-variable-single-component processes reuse the existing DOF
        // table.
        const bool manage_storage = false;
        return std::make_tuple(local_to_global_index_map_.get(),
                               manage_storage);
    }

    // Otherwise construct a new DOF table.
    std::vector<MeshLib::MeshSubset> all_mesh_subsets_single_component{
        *mesh_subset_all_nodes_};

    const bool manage_storage = true;
    return std::make_tuple(new NumLib::LocalToGlobalIndexMap(
                               std::move(all_mesh_subsets_single_component),
                               // by location order is needed for output
                               NumLib::ComponentOrder::BY_LOCATION),
                           manage_storage);
}

void HTProcess::setCoupledSolutionsOfPreviousTimeStepPerProcess(
    const int process_id)
{
    const auto& x_t0 = xs_previous_timestep_[process_id];
    if (x_t0 == nullptr)
    {
        OGS_FATAL(
            "Memory is not allocated for the global vector of the solution of "
            "the previous time step for the staggered scheme.\n It can be done "
            "by overriding Process::preTimestepConcreteProcess (ref. "
            "HTProcess::preTimestepConcreteProcess) ");
    }

    coupled_solutions_->coupled_xs_t0[process_id] = x_t0.get();
}

void HTProcess::setCoupledSolutionsOfPreviousTimeStep()
{
    coupled_solutions_->coupled_xs_t0.resize(2);
    setCoupledSolutionsOfPreviousTimeStepPerProcess(heat_transport_process_id_);
    setCoupledSolutionsOfPreviousTimeStepPerProcess(hydraulic_process_id_);
}

Eigen::Vector3d HTProcess::getFlux(std::size_t element_id,
                                   MathLib::Point3d const& p,
                                   double const t,
                                   std::vector<GlobalVector*> const& x) const
{
    // fetch local_x from primary variable
    std::vector<GlobalIndexType> indices_cache;
    auto const r_c_indices = NumLib::getRowColumnIndices(
        element_id, *local_to_global_index_map_, indices_cache);
    std::vector<std::vector<GlobalIndexType>> indices_of_all_coupled_processes{
        x.size(), r_c_indices.rows};
    auto const local_x =
        getCoupledLocalSolutions(x, indices_of_all_coupled_processes);

    return local_assemblers_[element_id]->getFlux(p, t, local_x);
}

// this is almost a copy of the implementation in the GroundwaterFlow
void HTProcess::postTimestepConcreteProcess(std::vector<GlobalVector*> const& x,
                                            const double t,
                                            const double /*delta_t*/,
                                            int const process_id)
{
    // For the monolithic scheme, process_id is always zero.
    if (use_monolithic_scheme_ && process_id != 0)
    {
        OGS_FATAL(
            "The condition of process_id = 0 must be satisfied for "
            "monolithic HTProcess, which is a single process.");
    }
    if (!use_monolithic_scheme_ && process_id != hydraulic_process_id_)
    {
        DBUG("This is the thermal part of the staggered HTProcess.");
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
}  // namespace HT
}  // namespace ProcessLib
