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

#include "NumLib/Extrapolation/LocalLinearLeastSquaresExtrapolator.h"
#include "ProcessLib/Process.h"
#include "GroundwaterFlowFEM.h"
#include "GroundwaterFlowProcessData.h"


namespace ProcessLib
{
namespace GroundwaterFlow
{

class GroundwaterFlowProcess final
        : public Process
{
    using Base = Process;

public:
    GroundwaterFlowProcess(
        MeshLib::Mesh& mesh,
        Base::NonlinearSolver& nonlinear_solver,
        std::unique_ptr<Base::TimeDiscretization>&& time_discretization,
        std::vector<std::reference_wrapper<ProcessVariable>>&& process_variables,
        GroundwaterFlowProcessData&& process_data,
        SecondaryVariableCollection&& secondary_variables,
        ProcessOutput&& process_output
        );

    //! \name ODESystem interface
    //! @{

    bool isLinear() const override
    {
        return true;
    }

    //! @}

private:
    using LocalAssemblerInterface =
        GroundwaterFlowLocalAssemblerInterface;

    using GlobalAssembler = NumLib::VectorMatrixAssembler<
            LocalAssemblerInterface,
            NumLib::ODESystemTag::FirstOrderImplicitQuasilinear>;

    using ExtrapolatorInterface = NumLib::Extrapolator<
        IntegrationPointValue, LocalAssemblerInterface>;
    using ExtrapolatorImplementation = NumLib::LocalLinearLeastSquaresExtrapolator<
        IntegrationPointValue, LocalAssemblerInterface>;

    void initializeConcreteProcess(
            NumLib::LocalToGlobalIndexMap const& dof_table,
            MeshLib::Mesh const& mesh,
            unsigned const integration_order) override;

    void assembleConcreteProcess(const double t, GlobalVector const& x,
                                 GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b) override;


    GroundwaterFlowProcessData _process_data;

    std::unique_ptr<GlobalAssembler> _global_assembler;
    std::vector<std::unique_ptr<LocalAssemblerInterface>> _local_assemblers;

    std::unique_ptr<ExtrapolatorInterface> _extrapolator;
};

std::unique_ptr<GroundwaterFlowProcess>
createGroundwaterFlowProcess(
    MeshLib::Mesh& mesh,
    Process::NonlinearSolver& nonlinear_solver,
    std::unique_ptr<Process::TimeDiscretization>&& time_discretization,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    BaseLib::ConfigTree const& config);

}   // namespace GroundwaterFlow
}   // namespace ProcessLib

#endif  // PROCESS_LIB_GROUNDWATERFLOWPROCESS_H_
