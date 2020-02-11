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

#include "CreateLiquidFlowMaterialProperties.h"
#include "LiquidFlowLocalAssembler.h"
#include "LiquidFlowMaterialProperties.h"
#include "MeshLib/PropertyVector.h"
// TODO(TF) used for output of flux, if output classes are ready this has to be changed
#include "MeshLib/IO/writeMeshToFile.h"
#include "ProcessLib/Utils/CreateLocalAssemblers.h"

namespace ProcessLib
{
namespace LiquidFlow
{
void checkMPLProperties(
    MeshLib::Mesh const& mesh,
    MaterialPropertyLib::MaterialSpatialDistributionMap const& media_map)
{
    DBUG("Check the media properties of LiquidFlow process ...");

    std::array const requiredPropertyMedium = {
        MaterialPropertyLib::PropertyType::porosity,
        MaterialPropertyLib::PropertyType::permeability};

    std::array const requiredPropertyLiquidPhase = {
        MaterialPropertyLib::PropertyType::viscosity,
        MaterialPropertyLib::PropertyType::density};

    std::array const requiredPropertySolidPhase = {
        MaterialPropertyLib::PropertyType::storage};

    for (auto const& element : mesh.getElements())
    {
        auto const element_id = element->getID();

        auto const& medium = *media_map.getMedium(element_id);
        MaterialPropertyLib::checkRequiredProperties(
            medium, requiredPropertyMedium);

        MaterialPropertyLib::checkRequiredProperties(
            medium.phase("AqueousLiquid"), requiredPropertyLiquidPhase);

        MaterialPropertyLib::checkRequiredProperties(
            medium.phase("Solid"), requiredPropertySolidPhase);
    }
    DBUG("Media properties verified.");
}

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
    MeshLib::PropertyVector<int> const* const material_ids,
    int const gravitational_axis_id,
    double const gravitational_acceleration,
    double const reference_temperature,
    BaseLib::ConfigTree const& config,
    std::unique_ptr<ProcessLib::SurfaceFluxData>&& surfaceflux)
    : Process(std::move(name), mesh, std::move(jacobian_assembler), parameters,
              integration_order, std::move(process_variables),
              std::move(secondary_variables)),
      _gravitational_axis_id(gravitational_axis_id),
      _gravitational_acceleration(gravitational_acceleration),
      _reference_temperature(reference_temperature),
      _material_properties(
          createLiquidFlowMaterialProperties(config, parameters, material_ids)),
      _process_data(std::move(process_data)),
      _surfaceflux(std::move(surfaceflux))
{
    DBUG("Create Liquid flow process.");
}

void LiquidFlowProcess::initializeConcreteProcess(
    NumLib::LocalToGlobalIndexMap const& dof_table,
    MeshLib::Mesh const& mesh,
    unsigned const integration_order)
{
    checkMPLProperties(mesh, *_process_data.media_map.get());

    const int process_id = 0;
    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];
    ProcessLib::createLocalAssemblers<LiquidFlowLocalAssembler>(
        mesh.getDimension(), mesh.getElements(), dof_table,
        pv.getShapeFunctionOrder(), _local_assemblers,
        mesh.isAxiallySymmetric(), integration_order, _gravitational_axis_id,
        _gravitational_acceleration, _reference_temperature,
        *_material_properties, _process_data);

    _secondary_variables.addSecondaryVariable(
        "darcy_velocity",
        makeExtrapolator(
            mesh.getDimension(), getExtrapolator(), _local_assemblers,
            &LiquidFlowLocalAssemblerInterface::getIntPtDarcyVelocity));
}

void LiquidFlowProcess::assembleConcreteProcess(
    const double t, double const dt, std::vector<GlobalVector*> const& x,
    std::vector<GlobalVector*> const& xdot, int const process_id,
    GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b)
{
    DBUG("Assemble LiquidFlowProcess.");

    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_table = {std::ref(*_local_to_global_index_map)};

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    // Call global assembler for each local assembly item.
    GlobalExecutor::executeSelectedMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assemble, _local_assemblers,
        pv.getActiveElementIDs(), dof_table, t, dt, x, xdot, process_id, M, K,
        b, _coupled_solutions);
}

void LiquidFlowProcess::assembleWithJacobianConcreteProcess(
    const double t, double const dt, std::vector<GlobalVector*> const& x,
    GlobalVector const& xdot, const double dxdot_dx, const double dx_dx,
    int const process_id, GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b,
    GlobalMatrix& Jac)
{
    DBUG("AssembleWithJacobian LiquidFlowProcess.");

    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_table = {std::ref(*_local_to_global_index_map)};
    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    // Call global assembler for each local assembly item.
    GlobalExecutor::executeSelectedMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assembleWithJacobian,
        _local_assemblers, pv.getActiveElementIDs(), dof_table, t, dt, x, xdot,
        dxdot_dx, dx_dx, process_id, M, K, b, Jac, _coupled_solutions);
}

void LiquidFlowProcess::computeSecondaryVariableConcrete(const double t,
                                                         GlobalVector const& x,
                                                         int const process_id)
{
    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    DBUG("Compute the velocity for LiquidFlowProcess.");
    GlobalExecutor::executeSelectedMemberOnDereferenced(
        &LiquidFlowLocalAssemblerInterface::computeSecondaryVariable,
        _local_assemblers, pv.getActiveElementIDs(), getDOFTable(process_id), t,
        x, _coupled_solutions);
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
        element_id, *_local_to_global_index_map, indices_cache);
    constexpr int process_id = 0;  // monolithic scheme.
    std::vector<double> local_x(x[process_id]->get(r_c_indices.rows));

    return _local_assemblers[element_id]->getFlux(p, t, local_x);
}

// this is almost a copy of the implementation in the GroundwaterFlow
void LiquidFlowProcess::postTimestepConcreteProcess(
    std::vector<GlobalVector*> const& x,
    const double t,
    const double /*dt*/,
    int const process_id)
{
    if (!_surfaceflux)  // computing the surfaceflux is optional
    {
        return;
    }

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];
    _surfaceflux->integrate(x, t, *this, process_id, _integration_order, _mesh,
                            pv.getActiveElementIDs());
    _surfaceflux->save(t);
}

}  // namespace LiquidFlow
}  // namespace ProcessLib
