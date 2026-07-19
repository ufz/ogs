// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "BHE_1U.h"

#include <numbers>

#include "BHEParameterValidation.h"
#include "FlowAndTemperatureControl.h"
#include "ParameterLib/SpatialPosition.h"
#include "Physics.h"
#include "ThermoMechanicalFlowProperties.h"

namespace ProcessLib
{
namespace HeatTransportBHE
{
namespace BHE
{
BHE_1U::BHE_1U(BoreholeGeometry const& borehole,
               RefrigerantProperties const& refrigerant,
               GroutParameters const& grout,
               FlowAndTemperatureControl const& flowAndTemperatureControl,
               PipeConfigurationUType const& pipes,
               bool const use_python_bcs)
    : BHECommonUType{borehole, refrigerant,   grout, flowAndTemperatureControl,
                     pipes,    use_python_bcs}
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

std::array<double, BHE_1U::number_of_unknowns> BHE_1U::pipeHeatCapacities()
    const
{
    double const rho_r = refrigerant.density;
    double const specific_heat_capacity = refrigerant.specific_heat_capacity;
    double const rho_g = grout.rho_g;
    double const porosity_g = grout.porosity_g;
    double const heat_cap_g = grout.heat_cap_g;

    return {{/*i1*/ rho_r * specific_heat_capacity,
             /*o1*/ rho_r * specific_heat_capacity,
             /*g1*/ (1.0 - porosity_g) * rho_g * heat_cap_g,
             /*g2*/ (1.0 - porosity_g) * rho_g * heat_cap_g}};
}

std::array<double, BHE_1U::number_of_unknowns> BHE_1U::pipeHeatConductions(
    ParameterLib::SpatialPosition const& /*pos*/) const
{
    double const lambda_r = refrigerant.thermal_conductivity;
    double const rho_r = refrigerant.density;
    double const Cp_r = refrigerant.specific_heat_capacity;
    double const alpha_L = _pipes.longitudinal_dispersion_length;
    double const porosity_g = grout.porosity_g;
    double const lambda_g = grout.lambda_g;

    double const velocity_norm = std::abs(flow_velocity);

    // Here we calculate the laplace coefficients in the governing
    // equations of BHE. These governing equations can be found in
    // 1) Diersch (2013) FEFLOW book on page 952, M.120-122, or
    // 2) Diersch (2011) Comp & Geosci 37:1122-1135, Eq. 19-22.
    auto const pipe_conduction =
        lambda_r + rho_r * Cp_r * alpha_L * velocity_norm;
    auto const grout_conduction = (1.0 - porosity_g) * lambda_g;
    return {{pipe_conduction,     // i1
             pipe_conduction,     // o1
             grout_conduction,    // g1
             grout_conduction}};  // g2
}

std::array<Eigen::Vector3d, BHE_1U::number_of_unknowns>
BHE_1U::pipeAdvectionVectors(Eigen::Vector3d const& elem_direction) const
{
    double const& rho_r = refrigerant.density;
    double const& Cp_r = refrigerant.specific_heat_capacity;

    auto const legs = flowLegs();
    auto leg_adv = [&](double const v_signed) -> Eigen::Vector3d
    { return rho_r * Cp_r * v_signed * elem_direction; };

    return {{leg_adv(legs[0]),  // i1
             leg_adv(legs[1]),  // o1
             {0, 0, 0},         // g1
             {0, 0, 0}}};       // g2
}

/// Thermal resistances due to grout-soil exchange.
///
/// Check if constraints regarding negative thermal resistances are violated
/// apply correction procedure.
/// Section (1.5.5) in FEFLOW White Papers Vol V.
std::array<double, 3> thermalResistancesGroutSoil(double const chi,
                                                  double const R_ar,
                                                  double const R_g)
{
    double R_gs = computeRgs(chi, R_g);
    double R_gg = computeRgg(chi, R_gs, R_ar, R_g);
    double new_chi = chi;

    auto constraint = [&]()
    { return 1.0 / ((1.0 / R_gg) + (1.0 / (2.0 * R_gs))); };

    std::array<double, 3> const multiplier{chi * 2.0 / 3.0, chi * 1.0 / 3.0,
                                           0.0};
    for (double m_chi : multiplier)
    {
        if (constraint() >= 0)
        {
            break;
        }
        DBUG(
            "Warning! Correction procedure was applied due to negative thermal "
            "resistance! Chi = {:f}.\n",
            m_chi);

        R_gs = computeRgs(m_chi, R_g);
        R_gg = computeRgg(m_chi, R_gs, R_ar, R_g);
        new_chi = m_chi;
    }

    return {new_chi, R_gg, R_gs};
}

void BHE_1U::updateHeatTransferCoefficients(double const flow_rate)
{
    auto const tm_flow = calculateThermoMechanicalFlowPropertiesPipe(
        _pipes.inlet, borehole_geometry.length, refrigerant, flow_rate);

    flow_velocity = tm_flow.velocity;
    cached_nu_ = tm_flow.nusselt_number;
}

std::vector<double> BHE_1U::thermalResistances(
    ParameterLib::SpatialPosition const& pos) const
{
    return calcThermalResistances(cached_nu_, pos);
}

/// Nu is the Nusselt number.
std::vector<double> BHE_1U::calcThermalResistances(
    double const Nu, ParameterLib::SpatialPosition const& pos) const
{
    constexpr double pi = std::numbers::pi;

    double const lambda_r = refrigerant.thermal_conductivity;
    double const lambda_g = grout.lambda_g;
    // t=0.0: borehole properties are physically time-invariant; genuinely
    // time-varying parameter types are rejected in createPipe.
    double const lambda_p =
        sampleStrictPositive(_pipes.inlet.wall_thermal_conductivity, 0.0, pos,
                             "wall_thermal_conductivity");

    // thermal resistances due to advective flow of refrigerant in the _pipes
    // Eq. 36 in Diersch_2011_CG
    double const R_adv_i1 = 1.0 / (Nu * lambda_r * pi);
    double const R_adv_o1 = 1.0 / (Nu * lambda_r * pi);

    // thermal resistance due to thermal conductivity of the pipe wall material
    // Eq. 49
    double const inlet_outside_diameter = _pipes.inlet.outsideDiameter();
    double const R_con_a = pipeWallThermalResistance(_pipes.inlet, lambda_p);

    // the average outer diameter of the _pipes
    double const d0 = inlet_outside_diameter;
    double const D = sampleStrictPositive(borehole_geometry.diameter, 0.0, pos,
                                          "borehole_diameter");
    // Context prefix shared by the acosh-argument checks below.
    auto const geometry_context = [&](std::string_view const equation)
    {
        return fmt::format(
            "BHE geometry invalid at {} ({:s}): D={:g}, d0={:g}, "
            "distance={:g}",
            pos, equation, D, d0, _pipes.distance);
    };
    // Eq. 51: chi requires D > sqrt(2) * d0 so that log(D/(sqrt(2)*d0)) > 0.
    checkBoreholeVsPipeDiameter(D, std::sqrt(2.0) * d0, pos,
                                "BHE 1U chi formula (Eq. 51)");
    double const chi = std::log(std::sqrt(D * D + 2 * d0 * d0) / 2 / d0) /
                       std::log(D / std::sqrt(2) / d0);
    // Eq. 52
    // thermal resistances of the grout
    double const acosh_arg_R_g =
        (D * D + d0 * d0 - _pipes.distance * _pipes.distance) / (2 * D * d0);
    checkAcoshArg(acosh_arg_R_g, geometry_context("BHE 1U R_g (Eq. 52)"));
    double const R_g = std::acosh(acosh_arg_R_g) / (2 * pi * lambda_g) *
                       (1.601 - 0.888 * _pipes.distance / D);

    // thermal resistance due to inter-grout exchange
    double const acosh_arg_R_ar =
        (2.0 * _pipes.distance * _pipes.distance - d0 * d0) / d0 / d0;
    checkAcoshArg(acosh_arg_R_ar, geometry_context("BHE 1U R_ar"));
    double const R_ar = std::acosh(acosh_arg_R_ar) / (2.0 * pi * lambda_g);

    auto const [chi_new, R_gg, R_gs] =
        thermalResistancesGroutSoil(chi, R_ar, R_g);

    // thermal resistance due to the grout transition.
    double const R_con_b = chi_new * R_g;
    // Eq. 29 and 30
    double const R_fig = R_adv_i1 + R_con_a + R_con_b;
    double const R_fog = R_adv_o1 + R_con_a + R_con_b;

    return {R_fig, R_fog, R_gg, R_gs};
}

std::array<std::pair<std::size_t /*node_id*/, int /*component*/>, 2>
BHE_1U::getBHEInflowDirichletBCNodesAndComponents(
    std::size_t const top_node_id,
    std::size_t const /*bottom_node_id*/,
    int const in_component_id)
{
    return {std::make_pair(top_node_id, in_component_id),
            std::make_pair(top_node_id, in_component_id + 1)};
}

std::optional<
    std::array<std::pair<std::size_t /*node_id*/, int /*component*/>, 2>>
BHE_1U::getBHEBottomDirichletBCNodesAndComponents(
    std::size_t const bottom_node_id,
    int const in_component_id,
    int const out_component_id)
{
    return {{std::make_pair(bottom_node_id, in_component_id),
             std::make_pair(bottom_node_id, out_component_id)}};
}

std::array<double, BHE_1U::number_of_unknowns> BHE_1U::crossSectionAreas(
    ParameterLib::SpatialPosition const& pos) const
{
    double const D = sampleStrictPositive(borehole_geometry.diameter, 0.0, pos,
                                          "borehole_diameter");
    double const borehole_area = Pipe::circleArea(D);
    double const half_borehole_area = borehole_area / number_of_grout_zones;
    return {
        {_pipes.inlet.area(), _pipes.outlet.area(),
         checkedGroutArea(half_borehole_area, _pipes.inlet.outsideArea(), pos),
         checkedGroutArea(half_borehole_area, _pipes.outlet.outsideArea(),
                          pos)}};
}

double BHE_1U::updateFlowRateAndTemperature(double const T_out,
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
