/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "HTProcess.h"

#include <cassert>

#include "MonolithicHTFEM.h"
#include "NumLib/DOF/DOFTableUtil.h"
#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "ProcessLib/SurfaceFlux/SurfaceFluxData.h"
#include "ProcessLib/Utils/CreateLocalAssemblers.h"
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
    std::unique_ptr<ProcessLib::SurfaceFluxData>&& surfaceflux)
    : Process(std::move(name), mesh, std::move(jacobian_assembler), parameters,
              integration_order, std::move(process_variables),
              std::move(secondary_variables), use_monolithic_scheme),
      _process_data(std::move(process_data)),
      _surfaceflux(std::move(surfaceflux))
{
}

void HTProcess::initializeConcreteProcess(
    NumLib::LocalToGlobalIndexMap const& dof_table,
    MeshLib::Mesh const& mesh,
    unsigned const integration_order)
{
    if (_use_monolithic_scheme)
    {
        ProcessLib::createLocalAssemblers<MonolithicHTFEM>(
            mesh.getDimension(), mesh.getElements(), dof_table,
            _local_assemblers, mesh.isAxiallySymmetric(), integration_order,
            _process_data);
    }
    else
    {
        ProcessLib::createLocalAssemblers<StaggeredHTFEM>(
            mesh.getDimension(), mesh.getElements(), dof_table,
            _local_assemblers, mesh.isAxiallySymmetric(), integration_order,
            _process_data);
    }

    _secondary_variables.addSecondaryVariable(
        "darcy_velocity",
        makeExtrapolator(mesh.getDimension(), getExtrapolator(),
                         _local_assemblers,
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
    if (_use_monolithic_scheme)
    {
        DBUG("Assemble HTProcess.");
        dof_tables.emplace_back(*_local_to_global_index_map);
    }
    else
    {
        if (process_id == _process_data.heat_transport_process_id)
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
        dof_tables.emplace_back(*_local_to_global_index_map);
        dof_tables.emplace_back(*_local_to_global_index_map);
    }

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];
    // Call global assembler for each local assembly item.
    GlobalExecutor::executeSelectedMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assemble, _local_assemblers,
        pv.getActiveElementIDs(), dof_tables, t, dt, x, xdot, process_id, M, K,
        b);
}

void HTProcess::assembleWithJacobianConcreteProcess(
    const double t, double const dt, std::vector<GlobalVector*> const& x,
    std::vector<GlobalVector*> const& xdot, int const process_id,
    GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b, GlobalMatrix& Jac)
{
    DBUG("AssembleWithJacobian HTProcess.");

    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_tables;
    if (!_use_monolithic_scheme)
    {
        dof_tables.emplace_back(std::ref(*_local_to_global_index_map));
    }
    else
    {
        dof_tables.emplace_back(std::ref(*_local_to_global_index_map));
        dof_tables.emplace_back(std::ref(*_local_to_global_index_map));
    }

    // Call global assembler for each local assembly item.
    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];
    GlobalExecutor::executeSelectedMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assembleWithJacobian,
        _local_assemblers, pv.getActiveElementIDs(), dof_tables, t, dt, x, xdot,
        process_id, M, K, b, Jac);
}

void HTProcess::setCoupledTermForTheStaggeredSchemeToLocalAssemblers(
    int const process_id)
{
    DBUG("Set the coupled term for the staggered scheme to local assembers.");

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];
    GlobalExecutor::executeSelectedMemberOnDereferenced(
        &HTLocalAssemblerInterface::setStaggeredCoupledSolutions,
        _local_assemblers, pv.getActiveElementIDs(), _coupled_solutions);
}

std::tuple<NumLib::LocalToGlobalIndexMap*, bool>
HTProcess::getDOFTableForExtrapolatorData() const
{
    if (!_use_monolithic_scheme)
    {
        // For single-variable-single-component processes reuse the existing DOF
        // table.
        const bool manage_storage = false;
        return std::make_tuple(_local_to_global_index_map.get(),
                               manage_storage);
    }

    // Otherwise construct a new DOF table.
    std::vector<MeshLib::MeshSubset> all_mesh_subsets_single_component{
        *_mesh_subset_all_nodes};

    const bool manage_storage = true;
    return std::make_tuple(new NumLib::LocalToGlobalIndexMap(
                               std::move(all_mesh_subsets_single_component),
                               // by location order is needed for output
                               NumLib::ComponentOrder::BY_LOCATION),
                           manage_storage);
}

Eigen::Vector3d HTProcess::getFlux(std::size_t element_id,
                                   MathLib::Point3d const& p,
                                   double const t,
                                   std::vector<GlobalVector*> const& x) const
{
    // fetch local_x from primary variable
    std::vector<GlobalIndexType> indices_cache;
    auto const r_c_indices = NumLib::getRowColumnIndices(
        element_id, *_local_to_global_index_map, indices_cache);
    std::vector<std::vector<GlobalIndexType>> indices_of_all_coupled_processes{
        x.size(), r_c_indices.rows};
    auto const local_x =
        getCoupledLocalSolutions(x, indices_of_all_coupled_processes);

    return _local_assemblers[element_id]->getFlux(p, t, local_x);
}

// this is almost a copy of the implementation in the GroundwaterFlow
void HTProcess::postTimestepConcreteProcess(std::vector<GlobalVector*> const& x,
                                            const double t,
                                            const double /*delta_t*/,
                                            int const process_id)
{
    // For the monolithic scheme, process_id is always zero.
    if (_use_monolithic_scheme && process_id != 0)
    {
        OGS_FATAL(
            "The condition of process_id = 0 must be satisfied for monolithic "
            "HTProcess, which is a single process.");
    }
    if (!_use_monolithic_scheme &&
        process_id != _process_data.hydraulic_process_id)
    {
        DBUG("This is the thermal part of the staggered HTProcess.");
        return;
    }
    if (!_surfaceflux)  // computing the surfaceflux is optional
    {
        return;
    }

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    _surfaceflux->integrate(x, t, *this, process_id, _integration_order, _mesh,
                            pv.getActiveElementIDs());
}
}  // namespace HT
}  // namespace ProcessLib
