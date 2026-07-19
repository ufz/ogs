// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "BHE_1P.h"

#include <numbers>

#include "BHEParameterValidation.h"
#include "FlowAndTemperatureControl.h"
#include "ParameterLib/SpatialPosition.h"
#include "Physics.h"
#include "ThermalResistanceHelpers.h"
#include "ThermoMechanicalFlowProperties.h"

namespace ProcessLib
{
namespace HeatTransportBHE
{
namespace BHE
{
BHE_1P::BHE_1P(BoreholeGeometry const& borehole,
               RefrigerantProperties const& refrigerant,
               GroutParameters const& grout,
               FlowAndTemperatureControl const& flowAndTemperatureControl,
               PipeConfiguration1PType const& pipes,
               bool const use_python_bcs)
    : BHECommon{borehole, refrigerant, grout, flowAndTemperatureControl,
                use_python_bcs},
      _pipe(pipes)
{
    // Initialize thermal resistances.
    auto values = visit(
        [&](auto const& control)
        {
            return control(refrigerant.reference_temperature,
                           0. /* initial time */);
        },
        flowAndTemperatureControl);
    updateHeatTransferCoefficients(values.flow_rate);
}

std::array<double, BHE_1P::number_of_unknowns> BHE_1P::pipeHeatCapacities()
    const
{
    double const rho_r = refrigerant.density;
    double const specific_heat_capacity = refrigerant.specific_heat_capacity;
    double const rho_g = grout.rho_g;
    double const porosity_g = grout.porosity_g;
    double const heat_cap_g = grout.heat_cap_g;

    return {{
        /*pipe*/ rho_r * specific_heat_capacity,
        /*grout*/ (1.0 - porosity_g) * rho_g * heat_cap_g,
    }};
}

std::array<double, BHE_1P::number_of_unknowns> BHE_1P::pipeHeatConductions(
    ParameterLib::SpatialPosition const& /*pos*/) const
{
    double const lambda_r = refrigerant.thermal_conductivity;
    double const rho_r = refrigerant.density;
    double const Cp_r = refrigerant.specific_heat_capacity;
    double const alpha_L = _pipe.longitudinal_dispersion_length;
    double const porosity_g = grout.porosity_g;
    double const lambda_g = grout.lambda_g;

    double const velocity_norm = std::abs(flow_velocity_);

    return {{
        // pipe
        (lambda_r + rho_r * Cp_r * alpha_L * velocity_norm),
        // grout
        (1.0 - porosity_g) * lambda_g,
    }};
}

std::array<Eigen::Vector3d, BHE_1P::number_of_unknowns>
BHE_1P::pipeAdvectionVectors(Eigen::Vector3d const& elem_direction) const
{
    double const& rho_r = refrigerant.density;
    double const& Cp_r = refrigerant.specific_heat_capacity;

    double const velocity = flow_velocity_;
    Eigen::Vector3d adv_vector = rho_r * Cp_r * velocity * elem_direction;

    return {// pipe
            adv_vector,
            // grout
            {0, 0, 0}};
}

void BHE_1P::updateHeatTransferCoefficients(double const flow_rate)
{
    auto const tm_flow = calculateThermoMechanicalFlowPropertiesPipe(
        _pipe.single_pipe, borehole_geometry.length, refrigerant, flow_rate);

    flow_velocity_ = tm_flow.velocity;
    cached_nu_ = tm_flow.nusselt_number;
}

std::vector<double> BHE_1P::thermalResistances(
    ParameterLib::SpatialPosition const& pos) const
{
    return calcThermalResistances(cached_nu_, pos);
}

std::vector<double> BHE_1P::calcThermalResistances(
    double const Nu, ParameterLib::SpatialPosition const& pos) const
{
    constexpr double pi = std::numbers::pi;

    double const lambda_r = refrigerant.thermal_conductivity;
    double const lambda_g = grout.lambda_g;
    // t=0.0: borehole properties are physically time-invariant (fixed at
    // drilling/installation); genuinely time-varying parameter types are
    // rejected in createPipe / createBoreholeGeometry.
    double const lambda_p =
        sampleStrictPositive(_pipe.single_pipe.wall_thermal_conductivity, 0.0,
                             pos, "wall_thermal_conductivity");

    // thermal resistances due to advective flow of refrigerant in the pipe
    double const R_adv_i1 = 1.0 / (Nu * lambda_r * pi);

    // thermal resistance due to thermal conductivity of the pipe wall material
    double const pipe_outside_diameter = _pipe.single_pipe.outsideDiameter();
    double const R_con_a =
        pipeWallThermalResistance(_pipe.single_pipe, lambda_p);

    // thermal resistances of the grout
    double const D = sampleStrictPositive(borehole_geometry.diameter, 0.0, pos,
                                          "borehole_diameter");
    checkBoreholeVsPipeDiameter(D, pipe_outside_diameter, pos,
                                "BHE 1P grout resistance");

    double const chi = std::log(std::sqrt(D * D + pipe_outside_diameter *
                                                      pipe_outside_diameter) /
                                std::sqrt(2) / pipe_outside_diameter) /
                       std::log(D / pipe_outside_diameter);
    double const R_g =
        std::log(D / pipe_outside_diameter) / 2 / (pi * lambda_g);

    double const R_con_b = chi * R_g;

    // thermal resistances due to grout-soil exchange
    double const R_gs = computeRgs(chi, R_g);

    double const R_fg = R_adv_i1 + R_con_a + R_con_b;

    return {R_fg, R_gs};
}

std::array<std::pair<std::size_t, int>, 2>
BHE_1P::getBHEInflowDirichletBCNodesAndComponents(
    std::size_t const top_node_id,
    std::size_t const bottom_node_id,
    int const in_component_id)
{
    return {std::make_pair(top_node_id, in_component_id),
            std::make_pair(bottom_node_id, in_component_id)};
}

std::optional<std::array<std::pair<std::size_t, int>, 2>>
BHE_1P::getBHEBottomDirichletBCNodesAndComponents(
    std::size_t const /*bottom_node_id*/,
    int const /*in_component_id*/,
    int const /*out_component_id*/)
{
    return {};
}

std::array<double, BHE_1P::number_of_unknowns> BHE_1P::crossSectionAreas(
    ParameterLib::SpatialPosition const& pos) const
{
    double const D = sampleStrictPositive(borehole_geometry.diameter, 0.0, pos,
                                          "borehole_diameter");
    double const borehole_area = Pipe::circleArea(D);
    return {{_pipe.single_pipe.area(),
             checkedGroutArea(borehole_area, _pipe.single_pipe.outsideArea(),
                              pos)}};
}

double BHE_1P::updateFlowRateAndTemperature(double const T_out,
                                            double const current_time)
{
    auto values =
        visit([&](auto const& control) { return control(T_out, current_time); },
              flowAndTemperatureControl);
    updateHeatTransferCoefficients(values.flow_rate);
    return values.temperature;
}
}  // namespace BHE
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
