/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "BHECommonCoaxial.h"

#include "Physics.h"
#include "ThermalResistancesCoaxial.h"
#include "ThermoMechanicalFlowProperties.h"

namespace ProcessLib
{
namespace HeatTransportBHE
{
namespace BHE
{
std::array<double, BHECommonCoaxial::number_of_unknowns>
BHECommonCoaxial::pipeHeatCapacities() const
{
    double const& rho_r = refrigerant.density;
    double const& specific_heat_capacity = refrigerant.specific_heat_capacity;
    double const& rho_g = grout.rho_g;
    double const& porosity_g = grout.porosity_g;
    double const& heat_cap_g = grout.heat_cap_g;

    return {{/*i*/ rho_r * specific_heat_capacity,
             /*o*/ rho_r * specific_heat_capacity,
             /*g*/ (1.0 - porosity_g) * rho_g * heat_cap_g}};
}

double BHECommonCoaxial::updateFlowRateAndTemperature(double const T_out,
                                             double const current_time)
{
    auto values = apply_visitor(
        [&](auto const& control) { return control(T_out, current_time); },
        flowAndTemperatureControl);
    updateHeatTransferCoefficients(values.flow_rate);
    return values.temperature;
}
std::array<double, BHECommonCoaxial::number_of_unknowns>
BHECommonCoaxial::pipeHeatConductions() const
{
    double const& lambda_r = refrigerant.thermal_conductivity;
    double const& rho_r = refrigerant.density;
    double const& Cp_r = refrigerant.specific_heat_capacity;
    double const& alpha_L = _pipes.longitudinal_dispersion_length;
    double const& porosity_g = grout.porosity_g;
    double const& lambda_g = grout.lambda_g;

    auto v = velocities();
    // Here we calculate the laplace coefficients in the governing
    // equations of BHE. These governing equations can be found in
    // 1) Diersch (2013) FEFLOW book on page 952, M.120-122, or
    // 2) Diersch (2011) Comp & Geosci 37:1122-1135, Eq. 26-28, 23-25
    return {{// pipe i, Eq. 26 and Eq. 23
             (lambda_r + rho_r * Cp_r * alpha_L * std::abs(v[0])),
             // pipe o, Eq. 27 and Eq. 24
             (lambda_r + rho_r * Cp_r * alpha_L * std::abs(v[1])),
             // pipe g, Eq. 28 and Eq. 25
             (1.0 - porosity_g) * lambda_g}};
}

std::array<Eigen::Vector3d, BHECommonCoaxial::number_of_unknowns>
BHECommonCoaxial::pipeAdvectionVectors() const
{
    double const& rho_r = refrigerant.density;
    double const& Cp_r = refrigerant.specific_heat_capacity;
    auto v = velocities();

    return {{// pipe i, Eq. 26 and Eq. 23
             {0, 0, -rho_r * Cp_r * std::abs(v[0])},
             // pipe o, Eq. 27 and Eq. 24
             {0, 0, rho_r * Cp_r * std::abs(v[1])},
             // grout g, Eq. 28 and Eq. 25
             {0, 0, 0}}};
}

std::array<double, BHECommonCoaxial::number_of_unknowns>
BHECommonCoaxial::calcThermalResistances(double const Nu_inner_pipe,
                                         double const Nu_annulus_pipe)
{
    // thermal resistances due to advective flow of refrigerant in the pipes
    auto const R_advective = calculateAdvectiveThermalResistance(
        _pipes.inner_pipe, _pipes.outer_pipe, refrigerant, Nu_inner_pipe,
        Nu_annulus_pipe);

    // thermal resistance due to thermal conductivity of the pipe wall material
    auto const R_conductive = calculatePipeWallThermalResistance(
        _pipes.inner_pipe, _pipes.outer_pipe);

    // thermal resistance due to the grout transition and grout-soil exchange.
    auto const R = calculateGroutAndGroutSoilExchangeThermalResistance(
        _pipes.outer_pipe, grout, borehole_geometry.diameter);

    // thermal resistance due to grout-soil exchange
    double const R_gs = R.grout_soil;

    // Eq. 56 and 57
    double const R_ff = R_advective.inner_pipe_coaxial + R_advective.a_annulus +
                        R_conductive.inner_pipe_coaxial;
    double const R_fg =
        R_advective.b_annulus + R_conductive.annulus + R.conductive_b;

    return getThermalResistances(R_gs, R_ff, R_fg);
}

void BHECommonCoaxial::updateHeatTransferCoefficients(double const flow_rate)
{
    auto const tm_flow_properties_annulus =
        calculateThermoMechanicalFlowPropertiesAnnulus(_pipes.inner_pipe,
                                                       _pipes.outer_pipe,
                                                       borehole_geometry.length,
                                                       refrigerant,
                                                       flow_rate);

    _flow_velocity_annulus = tm_flow_properties_annulus.velocity;

    auto const tm_flow_properties = calculateThermoMechanicalFlowPropertiesPipe(
        _pipes.inner_pipe, borehole_geometry.length, refrigerant, flow_rate);

    _flow_velocity_inner = tm_flow_properties.velocity;

    _thermal_resistances =
        calcThermalResistances(tm_flow_properties.nusselt_number,
                               tm_flow_properties_annulus.nusselt_number);
}
}  // namespace BHE
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
