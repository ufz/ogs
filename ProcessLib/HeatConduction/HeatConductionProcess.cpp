/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "HeatConductionProcess.h"

#include <cassert>

#include "ProcessLib/Utils/CreateLocalAssemblers.h"

namespace ProcessLib
{
namespace HeatConduction
{
HeatConductionProcess::HeatConductionProcess(
    MeshLib::Mesh& mesh,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    std::vector<std::reference_wrapper<ProcessVariable>>&& process_variables,
    HeatConductionProcessData&& process_data,
    SecondaryVariableCollection&& secondary_variables,
    NumLib::NamedFunctionCaller&& named_function_caller)
    : Process(mesh, parameters, std::move(process_variables),
              std::move(secondary_variables), std::move(named_function_caller)),
      _process_data(std::move(process_data))
{
}

void HeatConductionProcess::initializeConcreteProcess(
    NumLib::LocalToGlobalIndexMap const& dof_table,
    MeshLib::Mesh const& mesh,
    unsigned const integration_order)
{
    ProcessLib::createLocalAssemblers<LocalAssemblerData>(
        mesh.getDimension(), mesh.getElements(), dof_table, integration_order,
        _local_assemblers, _process_data);

    _secondary_variables.addSecondaryVariable(
        "heat_flux_x", 1,
        makeExtrapolator(
            getExtrapolator(), _local_assemblers,
            &HeatConductionLocalAssemblerInterface::getIntPtHeatFluxX));

    if (mesh.getDimension() > 1)
    {
        _secondary_variables.addSecondaryVariable(
            "heat_flux_y", 1,
            makeExtrapolator(
                getExtrapolator(), _local_assemblers,
                &HeatConductionLocalAssemblerInterface::getIntPtHeatFluxY));
    }
    if (mesh.getDimension() > 2)
    {
        _secondary_variables.addSecondaryVariable(
            "heat_flux_z", 1,
            makeExtrapolator(
                getExtrapolator(), _local_assemblers,
                &HeatConductionLocalAssemblerInterface::getIntPtHeatFluxZ));
    }
}

void HeatConductionProcess::assembleConcreteProcess(const double t,
                                                    GlobalVector const& x,
                                                    GlobalMatrix& M,
                                                    GlobalMatrix& K,
                                                    GlobalVector& b)
{
    DBUG("Assemble HeatConductionProcess.");

    // Call global assembler for each local assembly item.
    GlobalExecutor::executeMemberOnDereferenced(
        &HeatConductionLocalAssemblerInterface::assemble, _local_assemblers,
        *_local_to_global_index_map, t, x, M, K, b);
}

}  // namespace HeatConduction
}  // namespace ProcessLib
