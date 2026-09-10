// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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
///
/// \attention The process is restricted to vertical wells. Gravity enters the
/// balance equations twice, and the two paths do not treat an inclined well
/// consistently:
/// - The specific body force is projected onto the well axis, so it carries
///   the factor \f$\cos\theta\f$ with the inclination \f$\theta\f$ measured
///   from the vertical.
/// - The drift-flux closure is not projected. Both its drift flux velocity,
///   MaterialPropertyLib::driftFluxVelocity(), and its profile parameter,
///   MaterialPropertyLib::driftFluxProfileParameter(), are correlations for
///   vertical flow and receive no inclination correction.
///
/// A simulation on an inclined mesh therefore combines an inclination-corrected
/// body force with an uncorrected buoyant slip and is incorrect: the axial
/// drift, and with it the void fraction and the phase holdup, are
/// overestimated. The mismatch is worst for a horizontal section, where the
/// body force vanishes while the drift flux velocity keeps its full vertical
/// value. Beyond roughly 70 degrees from the vertical the flow leaves the
/// bubbly and slug regimes the closure is derived for, so an inclination
/// correction alone would not make such a section right either.
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
                                          GlobalVector const& x_prev,
                                          int const process_id) override;

    void postTimestepConcreteProcess(std::vector<GlobalVector*> const& x,
                                     std::vector<GlobalVector*> const& x_prev,
                                     const double t,
                                     const double dt,
                                     int const process_id) override;
};

}  // namespace WellboreSimulator
}  // namespace ProcessLib
