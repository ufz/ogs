// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <Eigen/Core>
#include <cmath>
#include <optional>

#include "BHECommon.h"
#include "BaseLib/Error.h"
#include "FlowAndTemperatureControl.h"
#include "PipeConfigurationCoaxial.h"

namespace ParameterLib
{
class SpatialPosition;
}

namespace ProcessLib
{
namespace HeatTransportBHE
{
namespace BHE
{
class BHECommonCoaxial : public BHECommon
{
public:
    BHECommonCoaxial(BoreholeGeometry const& borehole,
                     RefrigerantProperties const& refrigerant,
                     GroutParameters const& grout,
                     FlowAndTemperatureControl const& flowAndTemperatureControl,
                     PipeConfigurationCoaxial const& pipes,
                     bool const use_python_bcs);

    static constexpr int number_of_unknowns = 3;
    static constexpr int number_of_grout_zones = 1;
    static constexpr int number_of_flow_legs = 2;

    /// Signed fluid velocity per flow leg. Positive = flow along
    /// +elem_direction; negative = against it.  Leg ordering matches the
    /// BHE unknowns: { unknown 0 (inflow channel), unknown 1 (outflow
    /// channel) }; the physical-channel mapping (inner pipe vs annulus)
    /// is set per subclass via assignVelocities().
    std::array<double, number_of_flow_legs> flowLegs() const
    {
        return {+std::abs(velocity_inner), -std::abs(velocity_annulus)};
    }

    double updateFlowRateAndTemperature(double T_out, double current_time);

    std::vector<double> calcThermalResistances(
        double const Nu_inner_pipe,
        double const Nu_annulus_pipe,
        ParameterLib::SpatialPosition const& pos) const;

    /// Return the full vector of thermal resistances for the element at
    /// `pos`, computed once using the cached Nusselt numbers. Intended to
    /// be called once per element by the assembler so that the resistance
    /// computation is not repeated per unknown and per integration point.
    std::vector<double> thermalResistances(
        ParameterLib::SpatialPosition const& pos) const;

    std::array<double, number_of_unknowns> pipeHeatCapacities() const;

    static constexpr std::pair<int, int> inflow_outflow_bc_component_ids[] = {
        {0, 1}};

    static std::array<std::pair<std::size_t /*node_id*/, int /*component*/>, 2>
    getBHEInflowDirichletBCNodesAndComponents(
        std::size_t const top_node_id,
        std::size_t const /*bottom_node_id*/,
        int const in_component_id);

    static std::optional<
        std::array<std::pair<std::size_t /*node_id*/, int /*component*/>, 2>>
    getBHEBottomDirichletBCNodesAndComponents(std::size_t const bottom_node_id,
                                              int const in_component_id,
                                              int const out_component_id);

    std::array<double, number_of_unknowns> pipeHeatConductions() const;

    std::array<Eigen::Vector3d, number_of_unknowns> pipeAdvectionVectors(
        Eigen::Vector3d const& elem_direction) const;

    void updateHeatTransferCoefficients(double const flow_rate);

protected:
    PipeConfigurationCoaxial const _pipes;

    double cross_section_area_inner_pipe, cross_section_area_annulus;

    /// Fluid velocities indexed by unknown: `velocity_inner` is the
    /// velocity carried on unknown 0 (the inflow channel), and
    /// `velocity_annulus` is the velocity carried on unknown 1 (the
    /// outflow channel). The mapping from physical channels (inner pipe /
    /// annulus) to unknowns differs between CXA and CXC and is set once per
    /// flow-rate update via assignVelocities(). flowLegs() returns the
    /// signed per-leg velocities derived from these magnitudes.
    double velocity_inner = 0.0;
    double velocity_annulus = 0.0;
    double cached_nu_inner = 0.0;
    double cached_nu_annulus = 0.0;

    virtual std::vector<double> getThermalResistances(
        double const& R_gs, double const& R_ff, double const& R_fg) const = 0;

private:
    /// Assigns velocities from the physical channel velocities.
    /// Subclasses encode which channel maps to which unknown.
    virtual void assignVelocities(double inner_vel, double annulus_vel) = 0;
};
}  // namespace BHE
}  // namespace HeatTransportBHE
}  // end of namespace ProcessLib
