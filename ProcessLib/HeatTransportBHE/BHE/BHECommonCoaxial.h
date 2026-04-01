// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <Eigen/Core>
#include <optional>

#include "BHECommon.h"
#include "BaseLib/Error.h"
#include "FlowAndTemperatureControl.h"
#include "PipeConfigurationCoaxial.h"

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

    double updateFlowRateAndTemperature(double T_out, double current_time);

    std::vector<double> calcThermalResistances(
        double const Nu_inner_pipe,
        double const Nu_annulus_pipe,
        int const section_index = 0) const;

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

    std::array<double, number_of_unknowns> pipeHeatConductions(
        int const section_index = 0) const;

    std::array<Eigen::Vector3d, number_of_unknowns> pipeAdvectionVectors(
        Eigen::Vector3d const& /*elem_direction*/,
        int const section_index = 0) const;

    void updateHeatTransferCoefficients(double const flow_rate);

protected:
    PipeConfigurationCoaxial const _pipes;

    double cross_section_area_inner_pipe, cross_section_area_annulus;

    /// Returns fluid velocities indexed by unknown: {v_for_unknown_0,
    /// v_for_unknown_1}. The mapping from physical channels (inner pipe /
    /// annulus) to unknowns differs between CXA and CXC and is set once per
    /// flow-rate update via assignVelocities().
    std::array<double, 2> velocities() const
    {
        return {_flow_velocities[0], _flow_velocities[1]};
    }

    virtual std::vector<double> getThermalResistances(
        double const& R_gs, double const& R_ff, double const& R_fg) const = 0;

private:
    /// Assigns _flow_velocities from the physical channel velocities.
    /// Subclasses encode which channel maps to which unknown.
    virtual void assignVelocities(double inner_vel, double annulus_vel) = 0;
};
}  // namespace BHE
}  // namespace HeatTransportBHE
}  // end of namespace ProcessLib
