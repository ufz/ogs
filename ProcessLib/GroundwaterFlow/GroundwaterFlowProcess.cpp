/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "GroundwaterFlowProcess.h"

#include <cassert>

#include "ProcessLib/Utils/CreateLocalAssemblers.h"

namespace ProcessLib
{
namespace GroundwaterFlow
{

void checkMPLProperties(MeshLib::Mesh const& mesh,
                        GroundwaterFlowProcessData const& process_data)
{
    DBUG("Check the media properties of HT process ...");

    std::array const requiredPropertyMedium = {
        MaterialPropertyLib::PropertyType::reference_temperature,
        MaterialPropertyLib::PropertyType::diffusion};

    for (auto const& element : mesh.getElements())
    {
        auto const element_id = element->getID();

        auto const& medium = *process_data.media_map->getMedium(element_id);
        MaterialPropertyLib::checkRequiredProperties(
            medium, requiredPropertyMedium);
    }
    DBUG("Media properties verified.");
}

GroundwaterFlowProcess::GroundwaterFlowProcess(
    std::string name,
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    unsigned const integration_order,
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
        process_variables,
    GroundwaterFlowProcessData&& process_data,
    SecondaryVariableCollection&& secondary_variables,
    std::unique_ptr<ProcessLib::SurfaceFluxData>&& surfaceflux)
    : Process(std::move(name), mesh, std::move(jacobian_assembler), parameters,
              integration_order, std::move(process_variables),
              std::move(secondary_variables)),
      _process_data(std::move(process_data)),
      _surfaceflux(std::move(surfaceflux))
{
}

void GroundwaterFlowProcess::initializeConcreteProcess(
    NumLib::LocalToGlobalIndexMap const& dof_table,
    MeshLib::Mesh const& mesh,
    unsigned const integration_order)
{
    checkMPLProperties(mesh, _process_data);
    const int process_id = 0;
    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];
    ProcessLib::createLocalAssemblers<LocalAssemblerData>(
        mesh.getDimension(), mesh.getElements(), dof_table,
        pv.getShapeFunctionOrder(), _local_assemblers,
        mesh.isAxiallySymmetric(), integration_order, _process_data);

    _secondary_variables.addSecondaryVariable(
        "darcy_velocity",
        makeExtrapolator(
            mesh.getDimension(), getExtrapolator(), _local_assemblers,
            &GroundwaterFlowLocalAssemblerInterface::getIntPtDarcyVelocity));
}

void GroundwaterFlowProcess::assembleConcreteProcess(
    const double t, double const dt, std::vector<GlobalVector*> const& x,
    std::vector<GlobalVector*> const& xdot, int const process_id,
    GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b)
{
    DBUG("Assemble GroundwaterFlowProcess.");

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];
    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_table = {std::ref(*_local_to_global_index_map)};
    // Call global assembler for each local assembly item.
    GlobalExecutor::executeSelectedMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assemble, _local_assemblers,
        pv.getActiveElementIDs(), dof_table, t, dt, x, xdot, process_id, M, K,
        b, _coupled_solutions);
}

void GroundwaterFlowProcess::assembleWithJacobianConcreteProcess(
    const double t, double const dt, std::vector<GlobalVector*> const& x,
    GlobalVector const& xdot, const double dxdot_dx, const double dx_dx,
    int const process_id, GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b,
    GlobalMatrix& Jac)
{
    DBUG("AssembleWithJacobian GroundwaterFlowProcess.");

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];
    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_table = {std::ref(*_local_to_global_index_map)};
    // Call global assembler for each local assembly item.
    GlobalExecutor::executeSelectedMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assembleWithJacobian,
        _local_assemblers, pv.getActiveElementIDs(), dof_table, t, dt, x, xdot,
        dxdot_dx, dx_dx, process_id, M, K, b, Jac, _coupled_solutions);
}

}   // namespace GroundwaterFlow
}   // namespace ProcessLib
