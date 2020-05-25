/**
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file
 *
 * Created on August 19, 2016, 1:38 PM
 */

#include "LiquidFlowProcess.h"

#include <cassert>

#include "LiquidFlowLocalAssembler.h"
#include "MeshLib/PropertyVector.h"
// TODO(TF) used for output of flux, if output classes are ready this has to be changed
#include "MeshLib/IO/writeMeshToFile.h"
#include "ProcessLib/Utils/CreateLocalAssemblers.h"

namespace ProcessLib
{
namespace LiquidFlow
{
LiquidFlowProcess::LiquidFlowProcess(
    std::string name,
    MeshLib::Mesh& mesh,
    std::unique_ptr<AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    unsigned const integration_order,
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
        process_variables,
    LiquidFlowData&& process_data,
    SecondaryVariableCollection&& secondary_variables,
    std::unique_ptr<ProcessLib::SurfaceFluxData>&& surfaceflux)
    : Process(std::move(name), mesh, std::move(jacobian_assembler), parameters,
              integration_order, std::move(process_variables),
              std::move(secondary_variables)),
      process_data_(std::move(process_data)),
      surfaceflux_(std::move(surfaceflux))
{
    DBUG("Create Liquid flow process.");
}

void LiquidFlowProcess::initializeConcreteProcess(
    NumLib::LocalToGlobalIndexMap const& dof_table,
    MeshLib::Mesh const& mesh,
    unsigned const integration_order)
{
    const int process_id = 0;
    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];
    ProcessLib::createLocalAssemblers<LiquidFlowLocalAssembler>(
        mesh.getDimension(), mesh.getElements(), dof_table,
        pv.getShapeFunctionOrder(), local_assemblers_,
        mesh.isAxiallySymmetric(), integration_order, process_data_);

    secondary_variables_.addSecondaryVariable(
        "darcy_velocity",
        makeExtrapolator(
            mesh.getDimension(), getExtrapolator(), local_assemblers_,
            &LiquidFlowLocalAssemblerInterface::getIntPtDarcyVelocity));
}

void LiquidFlowProcess::assembleConcreteProcess(
    const double t, double const dt, std::vector<GlobalVector*> const& x,
    std::vector<GlobalVector*> const& xdot, int const process_id,
    GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b)
{
    DBUG("Assemble LiquidFlowProcess.");

    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_table = {std::ref(*local_to_global_index_map_)};

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    // Call global assembler for each local assembly item.
    GlobalExecutor::executeSelectedMemberDereferenced(
        global_assembler_, &VectorMatrixAssembler::assemble, local_assemblers_,
        pv.getActiveElementIDs(), dof_table, t, dt, x, xdot, process_id, M, K,
        b, coupled_solutions_);
}

void LiquidFlowProcess::assembleWithJacobianConcreteProcess(
    const double t, double const dt, std::vector<GlobalVector*> const& x,
    GlobalVector const& xdot, const double dxdot_dx, const double dx_dx,
    int const process_id, GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b,
    GlobalMatrix& Jac)
{
    DBUG("AssembleWithJacobian LiquidFlowProcess.");

    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_table = {std::ref(*local_to_global_index_map_)};
    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    // Call global assembler for each local assembly item.
    GlobalExecutor::executeSelectedMemberDereferenced(
        global_assembler_, &VectorMatrixAssembler::assembleWithJacobian,
        local_assemblers_, pv.getActiveElementIDs(), dof_table, t, dt, x, xdot,
        dxdot_dx, dx_dx, process_id, M, K, b, Jac, coupled_solutions_);
}

void LiquidFlowProcess::computeSecondaryVariableConcrete(
    double const t, double const dt, GlobalVector const& x,
    GlobalVector const& x_dot, int const process_id)
{
    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    DBUG("Compute the velocity for LiquidFlowProcess.");
    GlobalExecutor::executeSelectedMemberOnDereferenced(
        &LiquidFlowLocalAssemblerInterface::computeSecondaryVariable,
        local_assemblers_, pv.getActiveElementIDs(), getDOFTable(process_id), t,
        dt, x, x_dot, coupled_solutions_);
}

Eigen::Vector3d LiquidFlowProcess::getFlux(
    std::size_t const element_id,
    MathLib::Point3d const& p,
    double const t,
    std::vector<GlobalVector*> const& x) const
{
    // fetch local_x from primary variable
    std::vector<GlobalIndexType> indices_cache;
    auto const r_c_indices = NumLib::getRowColumnIndices(
        element_id, *local_to_global_index_map_, indices_cache);
    constexpr int process_id = 0;  // monolithic scheme.
    std::vector<double> local_x(x[process_id]->get(r_c_indices.rows));

    return local_assemblers_[element_id]->getFlux(p, t, local_x);
}

// this is almost a copy of the implementation in the GroundwaterFlow
void LiquidFlowProcess::postTimestepConcreteProcess(
    std::vector<GlobalVector*> const& x,
    const double t,
    const double /*dt*/,
    int const process_id)
{
    if (!surfaceflux_)  // computing the surfaceflux is optional
    {
        return;
    }

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];
    surfaceflux_->integrate(x, t, *this, process_id, integration_order_, mesh_,
                            pv.getActiveElementIDs());
    surfaceflux_->save(t);
}

}  // namespace LiquidFlow
}  // namespace ProcessLib
