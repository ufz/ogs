/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "GroundwaterFlowProcess.h"

#include <cassert>

#include "NumLib/Assembler/VectorMatrixAssembler.h"
#include "ProcessLib/Utils/CreateLocalAssemblers.h"

namespace ProcessLib
{
namespace GroundwaterFlow
{
GroundwaterFlowProcess::GroundwaterFlowProcess(
    MeshLib::Mesh& mesh,
    Base::NonlinearSolver& nonlinear_solver,
    std::unique_ptr<Base::TimeDiscretization>&& time_discretization,
    std::vector<std::reference_wrapper<ProcessVariable>>&& process_variables,
    GroundwaterFlowProcessData&& process_data,
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
            "GroundwaterFlowProcess can not be solved with the ForwardEuler"
            " time discretization scheme. Aborting");
        // Because the M matrix is not assembled. Thus, the linearized system
        // would be singular. The same applies to CrankNicolson with theta = 0.0,
        // but this case is not checked here.
        // Anyway, the GroundwaterFlowProcess shall be transferred to a simpler
        // ODESystemTag in the future.
    }
}

void GroundwaterFlowProcess::initializeConcreteProcess(
    NumLib::LocalToGlobalIndexMap const& dof_table,
    MeshLib::Mesh const& mesh,
    unsigned const integration_order)
{
    DBUG("Create global assembler.");
    _global_assembler.reset(new GlobalAssembler(dof_table));

    ProcessLib::createLocalAssemblers<LocalAssemblerData>(
        mesh.getDimension(), mesh.getElements(), dof_table, integration_order,
        _local_assemblers, _process_data);

    // TOOD Later on the DOF table can change during the simulation!
    _extrapolator.reset(new ExtrapolatorImplementation(
        Base::getMatrixSpecifications(), *Base::_local_to_global_index_map));

    Base::_secondary_variables.addSecondaryVariable(
        "darcy_velocity_x", 1,
        makeExtrapolator(IntegrationPointValue::DarcyVelocityX, *_extrapolator,
                         _local_assemblers));

    if (mesh.getDimension() > 1)
    {
        Base::_secondary_variables.addSecondaryVariable(
            "darcy_velocity_y", 1,
            makeExtrapolator(IntegrationPointValue::DarcyVelocityY,
                             *_extrapolator, _local_assemblers));
    }
    if (mesh.getDimension() > 2)
    {
        Base::_secondary_variables.addSecondaryVariable(
            "darcy_velocity_z", 1,
            makeExtrapolator(IntegrationPointValue::DarcyVelocityZ,
                             *_extrapolator, _local_assemblers));
    }
}

void GroundwaterFlowProcess::assembleConcreteProcess(const double t,
                                                     GlobalVector const& x,
                                                     GlobalMatrix& M,
                                                     GlobalMatrix& K,
                                                     GlobalVector& b)
{
    DBUG("Assemble GroundwaterFlowProcess.");

    // Call global assembler for each local assembly item.
    GlobalExecutor::executeMemberDereferenced(*_global_assembler,
                                              &GlobalAssembler::assemble,
                                              _local_assemblers, t, x, M, K, b);
}

}   // namespace GroundwaterFlow
}   // namespace ProcessLib
