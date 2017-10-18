/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file   LiquidFlowProcess.cpp
 *
 * Created on August 19, 2016, 1:38 PM
 */

#include "LiquidFlowProcess.h"

#include <cassert>

#include "MeshLib/PropertyVector.h"

#include "ProcessLib/Utils/CreateLocalAssemblers.h"
#include "LiquidFlowLocalAssembler.h"
#include "LiquidFlowMaterialProperties.h"

#include "MaterialLib/PorousMedium/Permeability/Permeability.h"
#include "MaterialLib/PorousMedium/Porosity/Porosity.h"
#include "MaterialLib/PorousMedium/Storage/Storage.h"

#include "CreateLiquidFlowMaterialProperties.h"

namespace ProcessLib
{
namespace LiquidFlow
{
LiquidFlowProcess::LiquidFlowProcess(
    MeshLib::Mesh& mesh,
    std::unique_ptr<AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order,
    std::vector<std::reference_wrapper<ProcessVariable>>&& process_variables,
    SecondaryVariableCollection&& secondary_variables,
    NumLib::NamedFunctionCaller&& named_function_caller,
    MeshLib::PropertyVector<int> const& material_ids,
    bool const has_material_ids,
    int const gravitational_axis_id,
    double const gravitational_acceleration,
    double const reference_temperature,
    BaseLib::ConfigTree const& config)
    : Process(mesh, std::move(jacobian_assembler), parameters,
              integration_order, std::move(process_variables),
              std::move(secondary_variables), std::move(named_function_caller)),
      _gravitational_axis_id(gravitational_axis_id),
      _gravitational_acceleration(gravitational_acceleration),
      _reference_temperature(reference_temperature),
      _material_properties(createLiquidFlowMaterialProperties(
          config, parameters, has_material_ids, material_ids))
{
    DBUG("Create Liquid flow process.");
}

void LiquidFlowProcess::initializeConcreteProcess(
    NumLib::LocalToGlobalIndexMap const& dof_table,
    MeshLib::Mesh const& mesh,
    unsigned const integration_order)
{
    ProcessLib::ProcessVariable const& pv = getProcessVariables()[0];
    ProcessLib::createLocalAssemblers<LiquidFlowLocalAssembler>(
        mesh.getDimension(), mesh.getElements(), dof_table,
        pv.getShapeFunctionOrder(), _local_assemblers,
        mesh.isAxiallySymmetric(), integration_order, _gravitational_axis_id,
        _gravitational_acceleration, _reference_temperature,
        *_material_properties);

    _secondary_variables.addSecondaryVariable(
        "darcy_velocity",
        makeExtrapolator(
            mesh.getDimension(), getExtrapolator(), _local_assemblers,
            &LiquidFlowLocalAssemblerInterface::getIntPtDarcyVelocity));
}

void LiquidFlowProcess::assembleConcreteProcess(const double t,
                                                GlobalVector const& x,
                                                GlobalMatrix& M,
                                                GlobalMatrix& K,
                                                GlobalVector& b)
{
    DBUG("Assemble LiquidFlowProcess.");
    // Call global assembler for each local assembly item.
    GlobalExecutor::executeMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assemble, _local_assemblers,
        *_local_to_global_index_map, t, x, M, K, b, _coupled_solutions);
}

void LiquidFlowProcess::assembleWithJacobianConcreteProcess(
    const double t, GlobalVector const& x, GlobalVector const& xdot,
    const double dxdot_dx, const double dx_dx, GlobalMatrix& M, GlobalMatrix& K,
    GlobalVector& b, GlobalMatrix& Jac)
{
    DBUG("AssembleWithJacobian LiquidFlowProcess.");

    // Call global assembler for each local assembly item.
    GlobalExecutor::executeMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assembleWithJacobian,
        _local_assemblers, *_local_to_global_index_map, t, x, xdot, dxdot_dx,
        dx_dx, M, K, b, Jac, _coupled_solutions);
}

void LiquidFlowProcess::computeSecondaryVariableConcrete(const double t,
                                                         GlobalVector const& x)
{
    DBUG("Compute the velocity for LiquidFlowProcess.");
    GlobalExecutor::executeMemberOnDereferenced(
        &LiquidFlowLocalAssemblerInterface::computeSecondaryVariable,
        _local_assemblers, *_local_to_global_index_map, t, x,
        _coupled_solutions);
}

}  // end of namespace
}  // end of namespace
