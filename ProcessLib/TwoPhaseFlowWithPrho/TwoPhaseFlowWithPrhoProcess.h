/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "ProcessLib/Process.h"
#include "TwoPhaseFlowWithPrhoLocalAssembler.h"

namespace MathLib
{
class PiecewiseLinearInterpolation;
}
namespace MeshLib
{
class Mesh;
}  // namespace MeshLib
namespace NumLib
{
class LocalToGlobalIndexMap;
}

namespace ProcessLib
{
namespace TwoPhaseFlowWithPrho
{
/**
 * \brief A class to simulate the two-phase flow process with P-rho model in
 * porous media
 */
class TwoPhaseFlowWithPrhoProcess final : public Process
{
public:
    TwoPhaseFlowWithPrhoProcess(
        std::string name,
        MeshLib::Mesh& mesh,
        std::unique_ptr<AbstractJacobianAssembler>&& jacobian_assembler,
        std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const&
            parameters,
        unsigned const integration_order,
        std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
            process_variables,
        TwoPhaseFlowWithPrhoProcessData&& process_data,
        SecondaryVariableCollection&& secondary_variables);

    bool isLinear() const override { return false; }

private:
    void initializeConcreteProcess(
        NumLib::LocalToGlobalIndexMap const& dof_table,
        MeshLib::Mesh const& mesh, unsigned const integration_order) override;

    void assembleConcreteProcess(const double t, double const dt,
                                 std::vector<GlobalVector*> const& x,
                                 std::vector<GlobalVector*> const& x_prev,
                                 int const process_id, GlobalMatrix& M,
                                 GlobalMatrix& K, GlobalVector& b) override;

    void assembleWithJacobianConcreteProcess(
        const double t, double const dt, std::vector<GlobalVector*> const& x,
        std::vector<GlobalVector*> const& x_prev, int const process_id,
        GlobalVector& b, GlobalMatrix& Jac) override;

    TwoPhaseFlowWithPrhoProcessData _process_data;

    std::vector<std::unique_ptr<TwoPhaseFlowWithPrhoLocalAssemblerInterface>>
        _local_assemblers;
};

}  // namespace TwoPhaseFlowWithPrho
}  // namespace ProcessLib
