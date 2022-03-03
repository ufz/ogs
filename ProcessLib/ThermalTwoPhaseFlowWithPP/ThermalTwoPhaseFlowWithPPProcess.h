/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"
#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "ProcessLib/Process.h"

#include "ThermalTwoPhaseFlowWithPPLocalAssembler.h"

namespace MeshLib
{
class Element;
class Mesh;
template <typename PROP_VAL_TYPE>
class PropertyVector;
}

namespace ProcessLib
{
namespace ThermalTwoPhaseFlowWithPP
{
/**
 * \brief A class to simulate the nonisothermal two-phase flow process with
 * phase change in porous media. Here the phase change is caused by
 * evaporation/condensation of the heavy component in the liquid phase.
 * Currently only the liquid phase is allowed to disappear due to the assumption
 * of a pure liquid phase. For the formulation of non-isothermal material
 * properties, see \cite Class_2002.
 */
class ThermalTwoPhaseFlowWithPPProcess final : public Process
{
public:
    ThermalTwoPhaseFlowWithPPProcess(
        std::string name,
        MeshLib::Mesh& mesh,
        std::unique_ptr<AbstractJacobianAssembler>&& jacobian_assembler,
        std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const&
            parameters,
        unsigned const integration_order,
        std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
            process_variables,
        ThermalTwoPhaseFlowWithPPProcessData&& process_data,
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

    void preTimestepConcreteProcess(std::vector<GlobalVector*> const& x,
                                    const double t, const double delta_t,
                                    const int process_id) override;

    ThermalTwoPhaseFlowWithPPProcessData _process_data;

    std::vector<
        std::unique_ptr<ThermalTwoPhaseFlowWithPPLocalAssemblerInterface>>
        _local_assemblers;
};

}  // namespace ThermalTwoPhaseFlowWithPP
}  // namespace ProcessLib
