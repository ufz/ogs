/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef OGS_TWOPHASEFLOWWITHPPPROCESS_H
#define OGS_TWOPHASEFLOWWITHPPPROCESS_H

#include "MaterialLib/TwoPhaseModels/TwoPhaseFlowWithPPMaterialProperties.h"
#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"
#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "ProcessLib/Process.h"
#include "TwoPhaseFlowWithPPLocalAssembler.h"

namespace MeshLib
{
class Element;
class Mesh;
template <typename PROP_VAL_TYPE>
class PropertyVector;
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
        MeshLib::Mesh& mesh,
        std::unique_ptr<AbstractJacobianAssembler>&& jacobian_assembler,
        std::vector<std::unique_ptr<ParameterBase>> const& parameters,
        unsigned const integration_order,
        std::vector<std::reference_wrapper<ProcessVariable>>&&
            process_variables,
        TwoPhaseFlowWithPPProcessData&& process_data,
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

    void assembleConcreteProcess(const double t, GlobalVector const& x,
                                 GlobalMatrix& M, GlobalMatrix& K,
                                 GlobalVector& b) override;

    void assembleWithJacobianConcreteProcess(
        const double t, GlobalVector const& x, GlobalVector const& xdot,
        const double dxdot_dx, const double dx_dx, GlobalMatrix& M,
        GlobalMatrix& K, GlobalVector& b, GlobalMatrix& Jac) override;

    TwoPhaseFlowWithPPProcessData _process_data;

    std::vector<std::unique_ptr<TwoPhaseFlowWithPPLocalAssemblerInterface>>
        _local_assemblers;
};

}  // end of namespace
}  // end of namespace

#endif /* TWOPHASEFLOWWITHPPPROCESS_H */
