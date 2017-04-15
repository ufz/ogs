/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "ProcessLib/Process.h"
#include "HeatConductionFEM.h"
#include "HeatConductionProcessData.h"

namespace ProcessLib
{
namespace HeatConduction
{
class HeatConductionProcess final : public Process
{

public:
    HeatConductionProcess(
        MeshLib::Mesh& mesh,
        std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&&
            jacobian_assembler,
        std::vector<std::unique_ptr<ParameterBase>> const& parameters,
        unsigned const integration_order,
        std::vector<std::reference_wrapper<ProcessVariable>>&&
            process_variables,
        HeatConductionProcessData&& process_data,
        SecondaryVariableCollection&& secondary_variables,
        NumLib::NamedFunctionCaller&& named_function_caller);

    //! \name ODESystem interface
    //! @{

    bool isLinear() const override { return true; }

    void computeSecondaryVariableConcrete(
        double const t, GlobalVector const& x,
        StaggeredCouplingTerm const& coupled_term) override;

    void preTimestepConcreteProcess(GlobalVector const& x, const double t,
                                    const double delta_t) override;

    // Get the solution of the previous time step.
    GlobalVector* getPreviousTimeStepSolution() const override
    {
        return _x_previous_timestep.get();
    }

private:
    void initializeConcreteProcess(
        NumLib::LocalToGlobalIndexMap const& dof_table,
        MeshLib::Mesh const& mesh,
        unsigned const integration_order) override;

    void assembleConcreteProcess(const double t, GlobalVector const& x,
                                 GlobalMatrix& M, GlobalMatrix& K,
                                 GlobalVector& b,
                                 StaggeredCouplingTerm const& coupling_term
                                ) override;

    void assembleWithJacobianConcreteProcess(
        const double t, GlobalVector const& x, GlobalVector const& xdot,
        const double dxdot_dx, const double dx_dx, GlobalMatrix& M,
        GlobalMatrix& K, GlobalVector& b, GlobalMatrix& Jac,
        StaggeredCouplingTerm const& coupling_term) override;

    HeatConductionProcessData _process_data;

    std::vector<std::unique_ptr<HeatConductionLocalAssemblerInterface>>
        _local_assemblers;

    /// Solution of the previous time step
    std::unique_ptr<GlobalVector> _x_previous_timestep = nullptr;
};

}  // namespace HeatConduction
}  // namespace ProcessLib
