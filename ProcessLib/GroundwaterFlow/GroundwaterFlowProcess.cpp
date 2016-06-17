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
#include "NumLib/Extrapolation/LocalLinearLeastSquaresExtrapolator.h"
#include "ProcessLib/Utils/CreateLocalAssemblers.h"

#include "GroundwaterFlowFEM.h"

namespace ProcessLib
{
namespace GroundwaterFlow
{
class GroundwaterFlowProcess::GroundwaterFlowProcessImpl
{
    using GlobalAssembler = NumLib::VectorMatrixAssembler<
        GroundwaterFlowLocalAssemblerInterface,
        NumLib::ODESystemTag::FirstOrderImplicitQuasilinear>;

    using ExtrapolatorInterface =
        NumLib::Extrapolator<IntegrationPointValue,
                             GroundwaterFlowLocalAssemblerInterface>;
    using ExtrapolatorImplementation =
        NumLib::LocalLinearLeastSquaresExtrapolator<
            IntegrationPointValue, GroundwaterFlowLocalAssemblerInterface>;

public:
    GroundwaterFlowProcessImpl(
        Process& process,
        NumLib::LocalToGlobalIndexMap& local_to_global_index_map,
        SecondaryVariableCollection& secondary_variables,
        GroundwaterFlowProcessData&& process_data)
        : _process(process),
          _local_to_global_index_map(local_to_global_index_map),
          _secondary_variables(secondary_variables),
          _process_data(std::move(process_data))
    {}


    void initializeConcreteProcess(
        NumLib::LocalToGlobalIndexMap const& dof_table,
        MeshLib::Mesh const& mesh,
        unsigned const integration_order)
    {
        DBUG("Create global assembler.");
        _global_assembler.reset(new GlobalAssembler(dof_table));

        ProcessLib::createLocalAssemblers<LocalAssemblerData>(
            mesh.getDimension(), mesh.getElements(), dof_table,
            integration_order, _local_assemblers, _process_data);

        // TOOD Later on the DOF table can change during the simulation!
        _extrapolator.reset(new ExtrapolatorImplementation(
            _process.getMatrixSpecifications(), _local_to_global_index_map));

        _secondary_variables.addSecondaryVariable(
            "darcy_velocity_x", 1,
            makeExtrapolator(IntegrationPointValue::DarcyVelocityX,
                             *_extrapolator, _local_assemblers));

        if (mesh.getDimension() > 1)
        {
            _secondary_variables.addSecondaryVariable(
                "darcy_velocity_y", 1,
                makeExtrapolator(IntegrationPointValue::DarcyVelocityY,
                                 *_extrapolator, _local_assemblers));
        }
        if (mesh.getDimension() > 2)
        {
            _secondary_variables.addSecondaryVariable(
                "darcy_velocity_z", 1,
                makeExtrapolator(IntegrationPointValue::DarcyVelocityZ,
                                 *_extrapolator, _local_assemblers));
        }
    }

    void assembleConcreteProcess(const double t, GlobalVector const& x,
                                 GlobalMatrix& M, GlobalMatrix& K,
                                 GlobalVector& b)
    {
        DBUG("Assemble GroundwaterFlowProcess.");

        // Call global assembler for each local assembly item.
        GlobalExecutor::executeMemberDereferenced(
            *_global_assembler, &GlobalAssembler::assemble, _local_assemblers,
            t, x, M, K, b);
    }

private:
    Process& _process;
    NumLib::LocalToGlobalIndexMap& _local_to_global_index_map;
    SecondaryVariableCollection& _secondary_variables;
    GroundwaterFlowProcessData _process_data;
    std::unique_ptr<GlobalAssembler> _global_assembler;
    std::vector<std::unique_ptr<GroundwaterFlowLocalAssemblerInterface>>
        _local_assemblers;
    std::unique_ptr<ExtrapolatorInterface> _extrapolator;
};

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
      impl(new GroundwaterFlowProcessImpl(
          *this, *Base::_local_to_global_index_map.get(),
          Base::_secondary_variables, std::move(process_data)))
{
    if (dynamic_cast<NumLib::ForwardEuler*>(
            &Base::getTimeDiscretization()) != nullptr)
    {
        OGS_FATAL(
            "GroundwaterFlowProcess can not be solved with the ForwardEuler"
            " time discretization scheme. Aborting");
        // Because the M matrix is not assembled. Thus, the linearized
        // system
        // would be singular. The same applies to CrankNicolson with theta =
        // 0.0,
        // but this case is not checked here.
        // Anyway, the GroundwaterFlowProcess shall be transferred to a
        // simpler
        // ODESystemTag in the future.
    }
}

// The destructor should be implemented after the implementation of
// GroundwaterFlowProcessImpl
GroundwaterFlowProcess::~GroundwaterFlowProcess()
{
}

void GroundwaterFlowProcess::initializeConcreteProcess(
    NumLib::LocalToGlobalIndexMap const& dof_table,
    MeshLib::Mesh const& mesh,
    unsigned const integration_order)
{
    impl->initializeConcreteProcess(dof_table, mesh, integration_order);
}

void GroundwaterFlowProcess::assembleConcreteProcess(const double t,
                                                     GlobalVector const& x,
                                                     GlobalMatrix& M,
                                                     GlobalMatrix& K,
                                                     GlobalVector& b)
{
    impl->assembleConcreteProcess(t, x, M, K, b);
}

std::unique_ptr<GroundwaterFlowProcess> createGroundwaterFlowProcess(
    MeshLib::Mesh& mesh,
    Process::NonlinearSolver& nonlinear_solver,
    std::unique_ptr<Process::TimeDiscretization>&& time_discretization,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{process__type}
    config.checkConfigParameter("type", "GROUNDWATER_FLOW");

    DBUG("Create GroundwaterFlowProcess.");

    // Process variable.
    auto process_variables = findProcessVariables(
        variables, config,
        {//! \ogs_file_param_special{process__GROUNDWATER_FLOW__process_variables__process_variable}
         "process_variable"});

    // Hydraulic conductivity parameter.
    auto& hydraulic_conductivity = findParameter<double,
                                                 MeshLib::Element const&>(
        config,
        //! \ogs_file_param_special{process__GROUNDWATER_FLOW__hydraulic_conductivity}
        "hydraulic_conductivity",
        parameters);

    DBUG("Use \'%s\' as hydraulic conductivity parameter.",
         hydraulic_conductivity.name.c_str());

    GroundwaterFlowProcessData process_data{hydraulic_conductivity};

    SecondaryVariableCollection secondary_variables{
        //! \ogs_file_param{process__secondary_variables}
        config.getConfigSubtreeOptional("secondary_variables"),
        {//! \ogs_file_param_special{process__GROUNDWATER_FLOW__secondary_variables__darcy_velocity_x}
         "darcy_velocity_x",
         //! \ogs_file_param_special{process__GROUNDWATER_FLOW__secondary_variables__darcy_velocity_y}
         "darcy_velocity_y",
         //! \ogs_file_param_special{process__GROUNDWATER_FLOW__secondary_variables__darcy_velocity_z}
         "darcy_velocity_z"}};

    ProcessOutput
        //! \ogs_file_param{process__output}
        process_output{config.getConfigSubtree("output"), process_variables,
                       secondary_variables};

    return std::unique_ptr<GroundwaterFlowProcess>{new GroundwaterFlowProcess{
        mesh, nonlinear_solver, std::move(time_discretization),
        std::move(process_variables), std::move(process_data),
        std::move(secondary_variables), std::move(process_output)}};
}

}  // namespace GroundwaterFlow
}  // namespace ProcessLib
