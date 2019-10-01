/**
 * \file
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "BHE_1P.h"

#include <boost/math/constants/constants.hpp>
#include "FlowAndTemperatureControl.h"
#include "Physics.h"
#include "ThermoMechanicalFlowProperties.h"

namespace ProcessLib
{
namespace HeatTransportBHE
{
namespace BHE
{
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

std::array<double, BHE_1P::number_of_unknowns> BHE_1P::pipeHeatConductions()
    const
{
    double const lambda_r = refrigerant.thermal_conductivity;
    double const rho_r = refrigerant.density;
    double const Cp_r = refrigerant.specific_heat_capacity;
    double const alpha_L = _pipes.longitudinal_dispersion_length;
    double const porosity_g = grout.porosity_g;
    double const lambda_g = grout.lambda_g;

    double const velocity_norm = std::abs(_flow_velocity) * std::sqrt(2);

    // Here we calculate the laplace coefficients in the governing
    // equations of BHE. These governing equations can be found in
    // 1) Diersch (2013) FEFLOW book on page 952, M.120-122, or
    // 2) Diersch (2011) Comp & Geosci 37:1122-1135, Eq. 19-22.
    return {{
        // pipe, Eq. 19
        (lambda_r + rho_r * Cp_r * alpha_L * velocity_norm),
        // grout, Eq. 21
        (1.0 - porosity_g) * lambda_g,
    }};
}

std::array<Eigen::Vector3d, BHE_1P::number_of_unknowns>
BHE_1P::pipeAdvectionVectors() const
{
    double const& rho_r = refrigerant.density;
    double const& Cp_r = refrigerant.specific_heat_capacity;

    return {{// pipe, Eq. 19
             {0, 0, -rho_r * Cp_r * _flow_velocity},
             // grout, Eq. 21
             {0, 0, 0}}};
}

double BHE_1P::compute_R_gs(double const chi, double const R_g)
{
    return (1 - chi) * R_g;
}

void BHE_1P::updateHeatTransferCoefficients(double const flow_rate)

{
    auto const tm_flow_properties = calculateThermoMechanicalFlowPropertiesPipe(
        _pipes.inlet, borehole_geometry.length, refrigerant, flow_rate);

    _flow_velocity = tm_flow_properties.velocity;
    _thermal_resistances =
        calcThermalResistances(tm_flow_properties.nusselt_number);
}

/// Nu is the Nusselt number.
std::array<double, BHE_1P::number_of_unknowns> BHE_1P::calcThermalResistances(
    double const Nu)
{
    constexpr double pi = boost::math::constants::pi<double>();

    double const lambda_r = refrigerant.thermal_conductivity;
    double const lambda_g = grout.lambda_g;
    double const lambda_p = _pipes.inlet.wall_thermal_conductivity;

    // thermal resistances due to advective flow of refrigerant in the _pipes
    // Eq. 36 in Diersch_2011_CG
    double const R_adv_i1 = 1.0 / (Nu * lambda_r * pi);

    // thermal resistance due to thermal conductivity of the pipe wall material
    // Eq. 49
    double const R_con_a =
        std::log(_pipes.outlet.diameter / _pipes.inlet.diameter) /
        (2.0 * pi * lambda_p);

    // the average outer diameter of the _pipes
    double const d0 = _pipes.outlet.diameter;
    double const D = borehole_geometry.diameter;
    // Eq. 51
    double const chi = std::log(std::sqrt(D * D + 2 * d0 * d0) / 2 / d0) /
                       std::log(D / std::sqrt(2) / d0);
    // Eq. 52
    // thermal resistances of the grout
    double const R_g =
        std::acosh((D * D + d0 * d0 - _pipes.distance * _pipes.distance) /
                   (2 * D * d0)) /
        (2 * pi * lambda_g) * (1.601 - 0.888 * _pipes.distance / D);

    // thermal resistance due to the grout transition.
    double const R_con_b = chi * R_g;
    // Eq. 29 and 30
    double const R_fig = R_adv_i1 + R_con_a + R_con_b;

    double R_gs;
    R_gs = compute_R_gs(chi, R_g);

    return {{R_fig, R_gs}};

    // keep the following lines------------------------------------------------
    // when debugging the code, printing the R and phi values are needed--------
    // std::cout << "Rfig =" << R_fig << " Rfog =" << R_fog << " Rgg =" <<
    // R_gg << " Rgs =" << R_gs << "\n"; double phi_fig = 1.0 / (R_fig *
    // S_i); double phi_fog = 1.0 / (R_fog * S_o); double phi_gg = 1.0 / (R_gg
    // * S_g1); double phi_gs = 1.0 / (R_gs * S_gs); std::cout << "phi_fig ="
    // << phi_fig << " phi_fog =" << phi_fog << " phi_gg =" << phi_gg << "
    // phi_gs =" << phi_gs << "\n";
    // -------------------------------------------------------------------------
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
