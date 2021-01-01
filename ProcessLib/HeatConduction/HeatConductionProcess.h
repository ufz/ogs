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
        std::string name,
        MeshLib::Mesh& mesh,
        std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&&
            jacobian_assembler,
        std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const&
            parameters,
        unsigned const integration_order,
        std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
            process_variables,
        HeatConductionProcessData&& process_data,
        SecondaryVariableCollection&& secondary_variables);

    //! \name ODESystem interface
    //! @{

    bool isLinear() const override { return true; }

    void computeSecondaryVariableConcrete(double const t, double const dt,
                                          std::vector<GlobalVector*> const& x,
                                          GlobalVector const& x_dot,
                                          int const process_id) override;

private:
    void initializeConcreteProcess(
        NumLib::LocalToGlobalIndexMap const& dof_table,
        MeshLib::Mesh const& mesh,
        unsigned const integration_order) override;

    void assembleConcreteProcess(const double t, double const /*dt*/,
                                 std::vector<GlobalVector*> const& x,
                                 std::vector<GlobalVector*> const& xdot,
                                 int const process_id, GlobalMatrix& M,
                                 GlobalMatrix& K, GlobalVector& b) override;

    void assembleWithJacobianConcreteProcess(
        const double t, double const /*dt*/,
        std::vector<GlobalVector*> const& x,
        std::vector<GlobalVector*> const& xdot, const double dxdot_dx,
        const double dx_dx, int const process_id, GlobalMatrix& M,
        GlobalMatrix& K, GlobalVector& b, GlobalMatrix& Jac) override;

    HeatConductionProcessData _process_data;

    std::vector<std::unique_ptr<HeatConductionLocalAssemblerInterface>>
        _local_assemblers;

    /// Solution of the previous time step
    std::unique_ptr<GlobalVector> _x_previous_timestep = nullptr;

    MeshLib::PropertyVector<double>* _heat_flux = nullptr;
};

}  // namespace HeatConduction
}  // namespace ProcessLib
