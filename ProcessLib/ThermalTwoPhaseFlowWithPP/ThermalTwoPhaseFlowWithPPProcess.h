/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
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
}  // namespace MeshLib

namespace ProcessLib
{
namespace ThermalTwoPhaseFlowWithPP
{
/**
 * A class to simulate the nonisothermal two-phase flow process in porous media
 * with max. 3 components. The water and air components are usually required on
 * a minimum. The air component is assumed to only exist in the gas phase. The
 * third component (organic contaminant) is optional, and can partition between
 * the liquid and gas phases. Technically, an algorithm for phase
 * appearance/disappearance is not implemented here. However, the liquid phase
 * is allowed to appear/disappear in some exceptional cases, e.g. the heat pipe
 * benchmark (since a pure liquid phase is assumed). For the formulation of
 * global conservation equations and constitutive relationships,
 * see \cite meng2021remediation.
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
