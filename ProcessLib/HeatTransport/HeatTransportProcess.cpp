
/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "HeatTransportProcess.h"

#include <cassert>

#include "ProcessLib/Utils/CreateLocalAssemblers.h"

namespace ProcessLib
{
namespace HeatTransport
{
HeatTransportProcess::HeatTransportProcess(
    MeshLib::Mesh& mesh,
    Base::NonlinearSolver& nonlinear_solver,
    std::unique_ptr<Base::TimeDiscretization>&& time_discretization,
    std::vector<std::reference_wrapper<ProcessVariable>>&& process_variables,
    HeatTransportProcessData&& process_data,
    SecondaryVariableCollection&& secondary_variables,
    ProcessOutput&& process_output)
    : Process(mesh, nonlinear_solver, std::move(time_discretization),
              std::move(process_variables), std::move(secondary_variables),
              std::move(process_output)),
      _process_data(std::move(process_data))
{
    if (dynamic_cast<NumLib::ForwardEuler*>(
            &Base::getTimeDiscretization()) != nullptr)
    {
        OGS_FATAL(
            "HeatTransportProcess can not be solved with the ForwardEuler"
            " time discretization scheme. Aborting");
        // Because the M matrix is not assembled. Thus, the linearized system
        // would be singular. The same applies to CrankNicolson with theta = 0.0,
        // but this case is not checked here.
        // Anyway, the HeatTransportProcess shall be transferred to a simpler
        // ODESystemTag in the future.
    }
}

void HeatTransportProcess::initializeConcreteProcess(
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
               &HeatTransportLocalAssemblerInterface::getIntPtHeatFluxX));

       if (mesh.getDimension() > 1) {
           _secondary_variables.addSecondaryVariable(
               "heat_flux_y", 1,
               makeExtrapolator(getExtrapolator(), _local_assemblers,
                                &HeatTransportLocalAssemblerInterface::
                                    getIntPtHeatFluxY));
       }
       if (mesh.getDimension() > 2) {
           _secondary_variables.addSecondaryVariable(
               "heat_flux_z", 1,
               makeExtrapolator(getExtrapolator(), _local_assemblers,
                                &HeatTransportLocalAssemblerInterface::
                                    getIntPtHeatFluxZ));
   }
}

void HeatTransportProcess::assembleConcreteProcess(const double t,
                                                     GlobalVector const& x,
                                                     GlobalMatrix& M,
                                                     GlobalMatrix& K,
                                                     GlobalVector& b)
{
    DBUG("Assemble HeatTransportProcess.");

    // Call global assembler for each local assembly item.
    GlobalExecutor::executeMemberOnDereferenced(
        &HeatTransportLocalAssemblerInterface::assemble,
        _local_assemblers, *_local_to_global_index_map, t, x, M, K, b);
}

}   // namespace HeatTransport
} // namespace ProcessLib
