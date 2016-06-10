/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESS_LIB_GROUNDWATERFLOWPROCESS_H_
#define PROCESS_LIB_GROUNDWATERFLOWPROCESS_H_

#include <cassert>

#include "NumLib/Assembler/VectorMatrixAssembler.h"
#include "NumLib/Extrapolation/LocalLinearLeastSquaresExtrapolator.h"
#include "ProcessLib/Process.h"
#include "ProcessLib/Utils/CreateLocalAssemblers.h"

#include "GroundwaterFlowFEM.h"
#include "GroundwaterFlowProcessData.h"


namespace ProcessLib
{
namespace GroundwaterFlow
{

template<typename GlobalSetup>
class GroundwaterFlowProcess final
        : public Process<GlobalSetup>
{
    using Base = Process<GlobalSetup>;
    using GlobalMatrix = typename GlobalSetup::MatrixType;
    using GlobalVector = typename GlobalSetup::VectorType;

public:
    GroundwaterFlowProcess(
        MeshLib::Mesh& mesh,
        typename Base::NonlinearSolver& nonlinear_solver,
        std::unique_ptr<typename Base::TimeDiscretization>&& time_discretization,
        std::vector<std::reference_wrapper<ProcessVariable>>&& process_variables,
        GroundwaterFlowProcessData&& process_data,
        SecondaryVariableCollection<GlobalVector>&& secondary_variables,
        ProcessOutput<GlobalVector>&& process_output
        )
        : Process<GlobalSetup>(mesh, nonlinear_solver, std::move(time_discretization),
                               std::move(process_variables),
                               std::move(secondary_variables),
                               std::move(process_output))
        , _process_data(std::move(process_data))
    {
        if (dynamic_cast<NumLib::ForwardEuler<GlobalVector>*>(
                    &Base::getTimeDiscretization()) != nullptr)
        {
            ERR("GroundwaterFlowProcess can not be solved with the ForwardEuler"
                " time discretization scheme. Aborting");
            // Because the M matrix is not assembled. Thus, the linearized system
            // would be singular. The same applies to CrankNicolson with theta = 0.0,
            // but this case is not checked here.
            // Anyway, the GroundwaterFlowProcess shall be transferred to a simpler
            // ODESystemTag in the future.
            std::abort();
        }
    }

    //! \name ODESystem interface
    //! @{

    bool isLinear() const override
    {
        return true;
    }

    //! @}

private:
    using LocalAssemblerInterface =
        GroundwaterFlowLocalAssemblerInterface<GlobalMatrix, GlobalVector>;

    using GlobalAssembler = NumLib::VectorMatrixAssembler<
            GlobalMatrix, GlobalVector, LocalAssemblerInterface,
            NumLib::ODESystemTag::FirstOrderImplicitQuasilinear>;

    using ExtrapolatorInterface = NumLib::Extrapolator<
        GlobalVector, IntegrationPointValue, LocalAssemblerInterface>;
    using ExtrapolatorImplementation = NumLib::LocalLinearLeastSquaresExtrapolator<
        GlobalVector, IntegrationPointValue, LocalAssemblerInterface>;

    void initializeConcreteProcess(
            NumLib::LocalToGlobalIndexMap const& dof_table,
            MeshLib::Mesh const& mesh,
            unsigned const integration_order) override
    {
        DBUG("Create global assembler.");
        _global_assembler.reset(new GlobalAssembler(dof_table));

        ProcessLib::createLocalAssemblers<GlobalSetup, LocalAssemblerData>(
                    mesh.getDimension(), mesh.getElements(),
                    dof_table, integration_order, _local_assemblers,
                    _process_data);

        // TOOD Later on the DOF table can change during the simulation!
        _extrapolator.reset(
            new ExtrapolatorImplementation(Base::getMatrixSpecifications(),
                                           *Base::_local_to_global_index_map));

        Base::_secondary_variables.addSecondaryVariable(
            "darcy_velocity_x", 1,
            makeExtrapolator(IntegrationPointValue::DarcyVelocityX, *_extrapolator,
                             _local_assemblers));

        if (mesh.getDimension() > 1) {
            Base::_secondary_variables.addSecondaryVariable(
                "darcy_velocity_y", 1,
                makeExtrapolator(IntegrationPointValue::DarcyVelocityY, *_extrapolator,
                                 _local_assemblers));
        }
        if (mesh.getDimension() > 2) {
            Base::_secondary_variables.addSecondaryVariable(
                "darcy_velocity_z", 1,
                makeExtrapolator(IntegrationPointValue::DarcyVelocityZ, *_extrapolator,
                                 _local_assemblers));
        }
    }

    void assembleConcreteProcess(const double t, GlobalVector const& x,
                                 GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b) override
    {
        DBUG("Assemble GroundwaterFlowProcess.");

        // Call global assembler for each local assembly item.
        GlobalSetup::executeMemberDereferenced(
            *_global_assembler, &GlobalAssembler::assemble,
            _local_assemblers, t, x, M, K, b);
    }


    GroundwaterFlowProcessData _process_data;

    std::unique_ptr<GlobalAssembler> _global_assembler;
    std::vector<std::unique_ptr<LocalAssemblerInterface>> _local_assemblers;

    std::unique_ptr<ExtrapolatorInterface> _extrapolator;
};

template <typename GlobalSetup>
std::unique_ptr<GroundwaterFlowProcess<GlobalSetup>>
createGroundwaterFlowProcess(
    MeshLib::Mesh& mesh,
    typename Process<GlobalSetup>::NonlinearSolver& nonlinear_solver,
    std::unique_ptr<typename Process<GlobalSetup>::TimeDiscretization>&& time_discretization,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{process__type}
    config.checkConfigParameter("type", "GROUNDWATER_FLOW");

    DBUG("Create GroundwaterFlowProcess.");

    // Process variable.
    auto process_variables = findProcessVariables(variables, config, {
        //! \ogs_file_param_special{process__GROUNDWATER_FLOW__process_variables__process_variable}
        "process_variable"
    });

    // Hydraulic conductivity parameter.
    auto& hydraulic_conductivity =
        findParameter<double, MeshLib::Element const&>(
            config,
            //! \ogs_file_param_special{process__GROUNDWATER_FLOW__hydraulic_conductivity}
            "hydraulic_conductivity",
            parameters);

    DBUG("Use \'%s\' as hydraulic conductivity parameter.",
         hydraulic_conductivity.name.c_str());

    GroundwaterFlowProcessData process_data {
        hydraulic_conductivity
    };

    SecondaryVariableCollection<typename GlobalSetup::VectorType> secondary_variables {
        //! \ogs_file_param{process__secondary_variables}
        config.getConfigSubtreeOptional("secondary_variables"),
        {
            //! \ogs_file_param_special{process__GROUNDWATER_FLOW__secondary_variables__darcy_velocity_x}
            "darcy_velocity_x",
            //! \ogs_file_param_special{process__GROUNDWATER_FLOW__secondary_variables__darcy_velocity_y}
            "darcy_velocity_y",
            //! \ogs_file_param_special{process__GROUNDWATER_FLOW__secondary_variables__darcy_velocity_z}
            "darcy_velocity_z"
        }
    };

    ProcessOutput<typename GlobalSetup::VectorType>
        //! \ogs_file_param{process__output}
        process_output{config.getConfigSubtree("output"),
                process_variables, secondary_variables};

    return std::unique_ptr<GroundwaterFlowProcess<GlobalSetup>>{
        new GroundwaterFlowProcess<GlobalSetup>{
            mesh, nonlinear_solver,std::move(time_discretization),
            std::move(process_variables),
            std::move(process_data),
            std::move(secondary_variables),
            std::move(process_output)
        }
    };
}

}   // namespace GroundwaterFlow
}   // namespace ProcessLib

#endif  // PROCESS_LIB_GROUNDWATERFLOWPROCESS_H_
