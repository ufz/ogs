/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "HTMaterialProperties.h"
#include "NumLib/Extrapolation/LocalLinearLeastSquaresExtrapolator.h"
#include "ProcessLib/Process.h"

namespace ProcessLib
{
namespace HT
{
class HTLocalAssemblerInterface;

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
              std::unique_ptr<HTMaterialProperties>&& material_properties,
              SecondaryVariableCollection&& secondary_variables,
              NumLib::NamedFunctionCaller&& named_function_caller);

    //! \name ODESystem interface
    //! @{

    bool isLinear() const override { return false; }
    //! @}

    // Get the solution of the previous time step.
    GlobalVector* getPreviousTimeStepSolution(
        const int process_id) const override
    {
        return _xs_previous_timestep[process_id].get();
    }

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

    void preTimestepConcreteProcess(
        GlobalVector const& x, double const t, double const dt,
        const int process_id) override;

    const std::unique_ptr<HTMaterialProperties> _material_properties;

    std::vector<std::unique_ptr<HTLocalAssemblerInterface>> _local_assemblers;
    
    /// Solutions of the previous time step
    std::array<std::unique_ptr<GlobalVector>, 2> _xs_previous_timestep;

};

}  // namespace HT
}  // namespace ProcessLib
