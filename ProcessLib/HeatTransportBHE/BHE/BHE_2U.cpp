// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "BHE_2U.h"

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
BHE_2U::BHE_2U(BoreholeGeometry const& borehole,
               RefrigerantProperties const& refrigerant,
               GroutParameters const& grout,
               FlowAndTemperatureControl const& flowAndTemperatureControl,
               PipeConfigurationUType const& pipes,
               bool const use_python_bcs)
    : BHECommonUType{borehole, refrigerant,   grout, flowAndTemperatureControl,
                     pipes,    use_python_bcs}
{
    checkEqualPipeOutsideDiameters(_pipes.inlet, _pipes.outlet, "BHE 2U");
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

std::array<double, BHE_2U::number_of_unknowns> BHE_2U::pipeHeatCapacities()
    const
{
    double const rho_r = refrigerant.density;
    double const specific_heat_capacity = refrigerant.specific_heat_capacity;
    double const rho_g = grout.rho_g;
    double const porosity_g = grout.porosity_g;
    double const heat_cap_g = grout.heat_cap_g;

    return {{/*i1*/ rho_r * specific_heat_capacity,
             /*i2*/ rho_r * specific_heat_capacity,
             /*o1*/ rho_r * specific_heat_capacity,
             /*o2*/ rho_r * specific_heat_capacity,
             /*g1*/ (1.0 - porosity_g) * rho_g * heat_cap_g,
             /*g2*/ (1.0 - porosity_g) * rho_g * heat_cap_g,
             /*g3*/ (1.0 - porosity_g) * rho_g * heat_cap_g,
             /*g4*/ (1.0 - porosity_g) * rho_g * heat_cap_g}};
}

std::array<double, BHE_2U::number_of_unknowns> BHE_2U::pipeHeatConductions()
    const
{
    double const lambda_r = refrigerant.thermal_conductivity;
    double const rho_r = refrigerant.density;
    double const Cp_r = refrigerant.specific_heat_capacity;
    double const alpha_L = _pipes.longitudinal_dispersion_length;
    double const porosity_g = grout.porosity_g;
    double const lambda_g = grout.lambda_g;

    double const velocity_norm = std::abs(flow_velocity);

    auto const pipe_conduction =
        lambda_r + rho_r * Cp_r * alpha_L * velocity_norm;
    auto const grout_conduction = (1.0 - porosity_g) * lambda_g;
    return {{pipe_conduction,     // i1
             pipe_conduction,     // i2
             pipe_conduction,     // o1
             pipe_conduction,     // o2
             grout_conduction,    // g1
             grout_conduction,    // g2
             grout_conduction,    // g3
             grout_conduction}};  // g4
}

std::array<Eigen::Vector3d, BHE_2U::number_of_unknowns>
BHE_2U::pipeAdvectionVectors(Eigen::Vector3d const& elem_direction) const
{
    double const rho_r = refrigerant.density;
    double const Cp_r = refrigerant.specific_heat_capacity;

    auto const legs = flowLegs();
    auto leg_adv = [&](double const v_signed) -> Eigen::Vector3d
    { return rho_r * Cp_r * v_signed * elem_direction; };

    return {{leg_adv(legs[0]),  // i1
             leg_adv(legs[1]),  // i2
             leg_adv(legs[2]),  // o1
             leg_adv(legs[3]),  // o2
             {0, 0, 0},         // g1
             {0, 0, 0},         // g2
             {0, 0, 0},         // g3
             {0, 0, 0}}};       // g4
}

/// Thermal resistances due to grout-soil exchange.
std::array<double, 4> thermalResistancesGroutSoil2U(double const chi,
                                                    double const R_ar_1,
                                                    double const R_ar_2,
                                                    double const R_g)
{
    double R_gs = computeRgs(chi, R_g);
    double R_gg_1 = computeRgg(chi, R_gs, R_ar_1, R_g);
    double R_gg_2 = computeRgg(chi, R_gs, R_ar_2, R_g);
    double chi_new = chi;

    auto constraint = [&]()
    { return 1.0 / ((1.0 / R_gg_1) + (1.0 / (2.0 * R_gs))); };

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
        R_gg_1 = computeRgg(m_chi, R_gs, R_ar_1, R_g);
        R_gg_2 = computeRgg(m_chi, R_gs, R_ar_2, R_g);
        chi_new = m_chi;
    }

    return {chi_new, R_gg_1, R_gg_2, R_gs};
}

void BHE_2U::updateHeatTransferCoefficients(double const flow_rate)
{
    auto const tm_flow = calculateThermoMechanicalFlowPropertiesPipe(
        _pipes.inlet, borehole_geometry.length, refrigerant, flow_rate);

    flow_velocity = tm_flow.velocity;
    cached_nu_ = tm_flow.nusselt_number;
}

std::vector<double> BHE_2U::thermalResistances(
    ParameterLib::SpatialPosition const& pos) const
{
    return calcThermalResistances(cached_nu_, pos);
}

std::vector<double> BHE_2U::calcThermalResistances(
    double const Nu, ParameterLib::SpatialPosition const& pos) const
{
    constexpr double pi = std::numbers::pi;

    double const lambda_r = refrigerant.thermal_conductivity;
    double const lambda_g = grout.lambda_g;
    // t=0.0: borehole properties are physically time-invariant; genuinely
    // time-varying parameter types are rejected in createPipe.
    double const lambda_p_inlet =
        sampleStrictPositive(_pipes.inlet.wall_thermal_conductivity, 0.0, pos,
                             "inlet wall_thermal_conductivity");
    double const lambda_p_outlet =
        sampleStrictPositive(_pipes.outlet.wall_thermal_conductivity, 0.0, pos,
                             "outlet wall_thermal_conductivity");

    // thermal resistances due to advective flow of refrigerant in the _pipes
    double const R_adv_i = 1.0 / (Nu * lambda_r * pi);
    double const R_adv_o = 1.0 / (Nu * lambda_r * pi);

    // thermal resistance due to thermal conductivity of the pipe wall material
    double const R_con_a_inlet =
        pipeWallThermalResistance(_pipes.inlet, lambda_p_inlet);
    double const R_con_a_outlet =
        pipeWallThermalResistance(_pipes.outlet, lambda_p_outlet);

    // Single pipe outside diameter used by the resistance formulas; the U-type
    // constructor enforces that inlet and outlet share it.
    double const d0 = _pipes.inlet.outsideDiameter();
    double const D = sampleStrictPositive(borehole_geometry.diameter, 0.0, pos,
                                          "borehole_diameter");
    // Context prefix shared by the acosh-argument checks below. `cause` names
    // the physical degeneracy that drives the argument to <= 1 so the fatal
    // message is actionable.
    auto const geometry_context =
        [&](std::string_view const equation, std::string_view const cause)
    {
        return uTypeGeometryContext(pos, D, d0, _pipes.distance, equation,
                                    cause);
    };
    // Eq. 38: chi requires D > 2 * d0 so that log(D/(2*d0)) > 0.
    checkBoreholeVsPipeDiameter(D, 2.0 * d0, pos,
                                "BHE 2U chi formula (Eq. 38)");
    double const chi =
        std::log(std::sqrt(D * D + 4 * d0 * d0) / 2 / std::sqrt(2) / d0) /
        std::log(D / 2 / d0);
    // Eq. 39
    double const acosh_arg_R_g =
        (D * D + d0 * d0 - 2 * _pipes.distance * _pipes.distance) /
        (2 * D * d0);
    checkAcoshArg(
        acosh_arg_R_g,
        geometry_context(
            "BHE 2U R_g (Eq. 39)",
            "The argument drops to <= 1 when a pipe reaches the borehole wall, "
            "giving a zero grout resistance R_g and an infinite (1/R) assembly "
            "coefficient."));
    double const R_g = std::acosh(acosh_arg_R_g) / (2 * pi * lambda_g) *
                       (3.098 - 4.432 * std::sqrt(2) * _pipes.distance / D +
                        2.364 * 2 * _pipes.distance * _pipes.distance / D / D);

    double const acosh_arg_R_ar_1 =
        (2.0 * _pipes.distance * _pipes.distance - d0 * d0) / d0 / d0;
    checkAcoshArg(
        acosh_arg_R_ar_1,
        geometry_context(
            "BHE 2U R_ar_1",
            "The argument drops to <= 1 when the pipes are too close, giving a "
            "zero inter-grout resistance R_ar and an infinite (1/R) assembly "
            "coefficient."));
    double const R_ar_1 = std::acosh(acosh_arg_R_ar_1) / (2.0 * pi * lambda_g);

    double const acosh_arg_R_ar_2 =
        (2.0 * 2.0 * _pipes.distance * _pipes.distance - d0 * d0) / d0 / d0;
    checkAcoshArg(
        acosh_arg_R_ar_2,
        geometry_context(
            "BHE 2U R_ar_2",
            "The argument drops to <= 1 when the pipes are too close, giving a "
            "zero inter-grout resistance R_ar and an infinite (1/R) assembly "
            "coefficient."));
    double const R_ar_2 = std::acosh(acosh_arg_R_ar_2) / (2.0 * pi * lambda_g);

    auto const [chi_new, R_gg_1, R_gg_2, R_gs] =
        thermalResistancesGroutSoil2U(chi, R_ar_1, R_ar_2, R_g);

    double const R_con_b = chi_new * R_g;

    double const R_fig = R_adv_i + R_con_a_inlet + R_con_b;
    double const R_fog = R_adv_o + R_con_a_outlet + R_con_b;

    return {R_fig, R_fog, R_gg_1, R_gg_2, R_gs};
}

std::array<std::pair<std::size_t, int>, 2>
BHE_2U::getBHEInflowDirichletBCNodesAndComponents(
    std::size_t const top_node_id,
    std::size_t const /*bottom_node_id*/,
    int const in_component_id)
{
    return {std::make_pair(top_node_id, in_component_id),
            std::make_pair(top_node_id, in_component_id + 2)};
}

std::optional<std::array<std::pair<std::size_t, int>, 2>>
BHE_2U::getBHEBottomDirichletBCNodesAndComponents(
    std::size_t const bottom_node_id,
    int const in_component_id,
    int const out_component_id)
{
    return {{std::make_pair(bottom_node_id, in_component_id),
             std::make_pair(bottom_node_id, out_component_id)}};
}

std::array<double, BHE_2U::number_of_unknowns> BHE_2U::crossSectionAreas(
    ParameterLib::SpatialPosition const& pos) const
{
    double const D = sampleStrictPositive(borehole_geometry.diameter, 0.0, pos,
                                          "borehole_diameter");
    double const borehole_area = Pipe::circleArea(D);
    double const quarter_borehole_area = borehole_area / number_of_grout_zones;
    double const grout_area_inlet = checkedGroutArea(
        quarter_borehole_area, _pipes.inlet.outsideArea(), pos);
    double const grout_area_outlet = checkedGroutArea(
        quarter_borehole_area, _pipes.outlet.outsideArea(), pos);

    return {{
        _pipes.inlet.area(),   // i1
        _pipes.inlet.area(),   // i2
        _pipes.outlet.area(),  // o1
        _pipes.outlet.area(),  // o2
        grout_area_inlet,      // g1
        grout_area_inlet,      // g2
        grout_area_outlet,     // g3
        grout_area_outlet,     // g4
    }};
}

double BHE_2U::updateFlowRateAndTemperature(double const T_out,
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
