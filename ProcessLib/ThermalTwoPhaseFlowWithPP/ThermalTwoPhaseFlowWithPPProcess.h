/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
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
 * \brief A class to simulate the nonisothermal two-phase flow process
 * with phase change phenomena between wetting phase and nonwetting phase in
 * porous media. Note that the phase change here can be caused by evaporation,
 * or by condensation, or by dissolution. The phase appearance and vanish
 * are also taken into account.
 */
class ThermalTwoPhaseFlowWithPPProcess final : public Process
{
public:
    ThermalTwoPhaseFlowWithPPProcess(
        MeshLib::Mesh& mesh,
        std::unique_ptr<AbstractJacobianAssembler>&& jacobian_assembler,
        std::vector<std::unique_ptr<ParameterBase>> const& parameters,
        unsigned const integration_order,
        std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
            process_variables,
        ThermalTwoPhaseFlowWithPPProcessData&& process_data,
        SecondaryVariableCollection&& secondary_variables,
        NumLib::NamedFunctionCaller&& named_function_caller,
        BaseLib::ConfigTree const& config,
        std::map<std::string,
                 std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
            curves);

    bool isLinear() const override { return false; }
private:
    void initializeConcreteProcess(
        NumLib::LocalToGlobalIndexMap const& dof_table,
        MeshLib::Mesh const& mesh, unsigned const integration_order) override;

    void assembleConcreteProcess(
        const double t, GlobalVector const& x, GlobalMatrix& M, GlobalMatrix& K,
        GlobalVector& b) override;

    void assembleWithJacobianConcreteProcess(
        const double t, GlobalVector const& x, GlobalVector const& xdot,
        const double dxdot_dx, const double dx_dx, GlobalMatrix& M,
        GlobalMatrix& K, GlobalVector& b, GlobalMatrix& Jac) override;

    void preTimestepConcreteProcess(GlobalVector const& x, const double t,
                                    const double delta_t,
                                    const int process_id) override;

    ThermalTwoPhaseFlowWithPPProcessData _process_data;

    std::vector<
        std::unique_ptr<ThermalTwoPhaseFlowWithPPLocalAssemblerInterface>>
        _local_assemblers;
};

}  // end of namespace
}  // end of namespace
