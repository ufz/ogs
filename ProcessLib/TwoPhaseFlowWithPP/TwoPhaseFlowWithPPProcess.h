/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "ProcessLib/Process.h"
#include "TwoPhaseFlowWithPPLocalAssembler.h"

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
namespace TwoPhaseFlowWithPP
{
/**
 * \brief A class to simulate the isothermal two-phase flow process with P-P
 * model in porous media.
 *
 * The gas and capillary pressure are used as primary variables.
 */
class TwoPhaseFlowWithPPProcess final : public Process
{
public:
    TwoPhaseFlowWithPPProcess(
        std::string name,
        MeshLib::Mesh& mesh,
        std::unique_ptr<AbstractJacobianAssembler>&& jacobian_assembler,
        std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const&
            parameters,
        unsigned const integration_order,
        std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
            process_variables,
        TwoPhaseFlowWithPPProcessData&& process_data,
        SecondaryVariableCollection&& secondary_variables,
        std::map<std::string,
                 std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
            curves);

    bool isLinear() const override { return false; }

private:
    void initializeConcreteProcess(
        NumLib::LocalToGlobalIndexMap const& dof_table,
        MeshLib::Mesh const& mesh, unsigned const integration_order) override;

    void assembleConcreteProcess(const double t, double const dt,
                                 std::vector<GlobalVector*> const& x,
                                 std::vector<GlobalVector*> const& xdot,
                                 int const process_id, GlobalMatrix& M,
                                 GlobalMatrix& K, GlobalVector& b) override;

    void assembleWithJacobianConcreteProcess(
        const double t, double const dt, std::vector<GlobalVector*> const& x,
        std::vector<GlobalVector*> const& xdot, int const process_id,
        GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b,
        GlobalMatrix& Jac) override;

    TwoPhaseFlowWithPPProcessData _process_data;

    std::vector<std::unique_ptr<TwoPhaseFlowWithPPLocalAssemblerInterface>>
        _local_assemblers;
};

}  // namespace TwoPhaseFlowWithPP
}  // namespace ProcessLib
