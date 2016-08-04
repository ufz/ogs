/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESS_LIB_HEATTRANSPORTPROCESS_H_
#define PROCESS_LIB_HEATTRANSPORTPROCESS_H_

#include "NumLib/Extrapolation/LocalLinearLeastSquaresExtrapolator.h"
#include "ProcessLib/Process.h"
#include "HeatTransportFEM.h"
#include "HeatTransportProcessData.h"

namespace ProcessLib
{
namespace HeatTransport
{
class HeatTransportProcess final : public Process
{
    using Base = Process;

public:
    HeatTransportProcess(
        MeshLib::Mesh& mesh,
        Base::NonlinearSolver& nonlinear_solver,
        std::unique_ptr<Base::TimeDiscretization>&& time_discretization,
        std::vector<std::reference_wrapper<ProcessVariable>>&&
            process_variables,
        HeatTransportProcessData&& process_data,
        SecondaryVariableCollection&& secondary_variables,
        ProcessOutput&& process_output);

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

    HeatTransportProcessData _process_data;

    std::vector<std::unique_ptr<HeatTransportLocalAssemblerInterface>>
        _local_assemblers;
};

}  // namespace HeatTransport
}  // namespace ProcessLib

#endif  // PROCESS_LIB_HEATTRANSPORTPROCESS_H_
