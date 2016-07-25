/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESS_LIB_HEATCONDUCTIONPROCESS_H_
#define PROCESS_LIB_HEATCONDUCTIONPROCESS_H_

#include "NumLib/Extrapolation/LocalLinearLeastSquaresExtrapolator.h"
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
        std::vector<std::unique_ptr<ParameterBase>> const& parameters,
        std::vector<std::reference_wrapper<ProcessVariable>>&&
            process_variables,
        HeatConductionProcessData&& process_data,
        SecondaryVariableCollection&& secondary_variables,
        NumLib::NamedFunctionCaller&& named_function_caller);

    //! \name ODESystem interface
    //! @{

    bool isLinear() const override { return true; }
    //! @}

private:
    void initializeConcreteProcess(
        NumLib::LocalToGlobalIndexMap const& dof_table,
        MeshLib::Mesh const& mesh,
        unsigned const integration_order) override;

    void assembleConcreteProcess(const double t, GlobalVector const& x,
                                 GlobalMatrix& M, GlobalMatrix& K,
                                 GlobalVector& b) override;

    HeatConductionProcessData _process_data;

    std::vector<std::unique_ptr<HeatConductionLocalAssemblerInterface>>
        _local_assemblers;
};

}  // namespace HeatConduction
}  // namespace ProcessLib

#endif  // PROCESS_LIB_HEATCONDUCTIONPROCESS_H_
