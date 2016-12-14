/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "HTFEM.h"
#include "HTProcessData.h"
#include "NumLib/Extrapolation/LocalLinearLeastSquaresExtrapolator.h"
#include "ProcessLib/Process.h"

namespace ProcessLib
{
namespace HT
{
/**
 * # HT process
 *
 * The implementation uses a monolithic approach, i.e., both processes
 * are assembled within one global system of equations.
 *
 * ## Process Coupling
 *
 * The advective term of the heat conduction equation is given by the confined
 * groundwater flow process, i.e., the heat conduction depends on darcy velocity
 * of the groundwater flow process. On the other hand the temperature
 * dependencies of the viscosity and density in the groundwater flow couples the
 * H process to the T process.
 *
 * \note At the moment there is not any coupling by source or sink terms, i.e.,
 * the coupling is implemented only by density changes due to temperature
 * changes in the buoyance term in the groundwater flow. This coupling schema is
 * referred to as the Boussinesq approximation.
 * */
class HTProcess final : public Process
{
public:
    HTProcess(MeshLib::Mesh& mesh,
              std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&&
                  jacobian_assembler,
              std::vector<std::unique_ptr<ParameterBase>> const& parameters,
              unsigned const integration_order,
              std::vector<std::reference_wrapper<ProcessVariable>>&&
                  process_variables,
              HTProcessData&& process_data,
              SecondaryVariableCollection&& secondary_variables,
              NumLib::NamedFunctionCaller&& named_function_caller);

    //! \name ODESystem interface
    //! @{

    bool isLinear() const override { return false; }
    //! @}

    ProcessType getProcessType() const override
                     {return ProcessLib::ProcessType::HTProcess;}

private:
    void initializeConcreteProcess(
        NumLib::LocalToGlobalIndexMap const& dof_table,
        MeshLib::Mesh const& mesh,
        unsigned const integration_order) override;

    void assembleConcreteProcess(const double t, GlobalVector const& x,
                                 GlobalMatrix& M, GlobalMatrix& K,
                                 GlobalVector& b) override;

    void assembleWithJacobianConcreteProcess(
        const double t, GlobalVector const& x, GlobalVector const& xdot,
        const double dxdot_dx, const double dx_dx, GlobalMatrix& M,
        GlobalMatrix& K, GlobalVector& b, GlobalMatrix& Jac) override;

    HTProcessData _process_data;

    std::vector<std::unique_ptr<HTLocalAssemblerInterface>> _local_assemblers;
};

}  // namespace HT
}  // namespace ProcessLib
