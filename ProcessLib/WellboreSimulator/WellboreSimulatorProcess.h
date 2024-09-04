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

#include "ProcessLib/Process.h"
#include "WellboreSimulatorProcessData.h"

namespace NumLib
{
class LocalToGlobalIndexMap;
}

namespace ProcessLib
{
namespace WellboreSimulator
{
class WellboreSimulatorLocalAssemblerInterface;

/// Under the assumptions of the drift-flux model in a geothermal well with
/// a constant cross-sectional area, the transient geo-fluid flow in a two-phase
/// one-component geothermal well can be quantified using one-dimensional
/// formulations of mass, momentum, and energy balance.
class WellboreSimulatorProcess final : public Process
{
public:
    WellboreSimulatorProcess(
        std::string name,
        MeshLib::Mesh& mesh,
        std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&&
            jacobian_assembler,
        std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const&
            parameters,
        unsigned const integration_order,
        std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
            process_variables,
        WellboreSimulatorProcessData&& process_data,
        SecondaryVariableCollection&& secondary_variables);
    //! \name ODESystem interface
    //! @{

    bool isLinear() const override { return false; }
    //! @}

private:
    void initializeConcreteProcess(
        NumLib::LocalToGlobalIndexMap const& dof_table,
        MeshLib::Mesh const& mesh,
        unsigned const integration_order) override;

    void assembleConcreteProcess(const double t, double const dt,
                                 std::vector<GlobalVector*> const& x,
                                 std::vector<GlobalVector*> const& x_prev,
                                 int const process_id, GlobalMatrix& M,
                                 GlobalMatrix& K, GlobalVector& b) override;

    void assembleWithJacobianConcreteProcess(
        const double t, double const dt, std::vector<GlobalVector*> const& x,
        std::vector<GlobalVector*> const& x_prev, int const process_id,
        GlobalVector& b, GlobalMatrix& Jac) override;

    WellboreSimulatorProcessData _process_data;

    std::vector<std::unique_ptr<WellboreSimulatorLocalAssemblerInterface>>
        _local_assemblers;

    void computeSecondaryVariableConcrete(double const t,
                                          double const dt,
                                          std::vector<GlobalVector*> const& x,
                                          GlobalVector const& x_dot,
                                          int const process_id) override;

    void postTimestepConcreteProcess(std::vector<GlobalVector*> const& x,
                                     std::vector<GlobalVector*> const& x_dot,
                                     const double t,
                                     const double dt,
                                     int const process_id) override;
};

}  // namespace WellboreSimulator
}  // namespace ProcessLib
