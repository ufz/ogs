/**
 * \file
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "BHE_2U.h"

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

    // Here we calculate the laplace coefficients in the governing
    // equations of BHE. These governing equations can be found in
    // 1) Diersch (2013) FEFLOW book on page 952, M.120-122, or
    // 2) Diersch (2011) Comp & Geosci 37:1122-1135, Eq. 19-22.
    return {{// pipe i1
             (lambda_r + rho_r * Cp_r * alpha_L * _flow_velocity),
             // pipe i2
             (lambda_r + rho_r * Cp_r * alpha_L * _flow_velocity),
             // pipe o1
             (lambda_r + rho_r * Cp_r * alpha_L * _flow_velocity),
             // pipe o2
             (lambda_r + rho_r * Cp_r * alpha_L * _flow_velocity),
             // pipe g1
             (1.0 - porosity_g) * lambda_g,
             // pipe g2
             (1.0 - porosity_g) * lambda_g,
             // pipe g3
             (1.0 - porosity_g) * lambda_g,
             // pipe g4
             (1.0 - porosity_g) * lambda_g}};
}

std::array<Eigen::Vector3d, BHE_2U::number_of_unknowns>
BHE_2U::pipeAdvectionVectors() const
{
    double const rho_r = refrigerant.density;
    double const Cp_r = refrigerant.specific_heat_capacity;

    return {{// pipe i1
             {0, 0, -rho_r * Cp_r * _flow_velocity},
             // pipe i2
             {0, 0, -rho_r * Cp_r * _flow_velocity},
             // pipe o1
             {0, 0, rho_r * Cp_r * _flow_velocity},
             // pipe o2
             {0, 0, rho_r * Cp_r * _flow_velocity},
             // grout g1
             {0, 0, 0},
             // grout g2
             {0, 0, 0},
             // grout g3
             {0, 0, 0},
             // grout g4
             {0, 0, 0}}};
}

double compute_R_gs_2U(double const chi, double const R_g)
{
    return (1 - chi) * R_g;
}

double compute_R_gg_2U(double const chi, double const R_gs, double const R_ar,
                       double const R_g)
{
    double const R_gg = 2.0 * R_gs * (R_ar - 2.0 * chi * R_g) /
                        (2.0 * R_gs - R_ar + 2.0 * chi * R_g);
    if (!std::isfinite(R_gg))
    {
        OGS_FATAL(
            "Error!!! Grout Thermal Resistance is an infinite number! The "
            "simulation will be stopped!");
    }

    return R_gg;
}

/// Thermal resistances due to grout-soil exchange.
///
/// Check if constraints regarding negative thermal resistances are violated
/// apply correction procedure.
/// Section (1.5.5) in FEFLOW White Papers Vol V.
std::vector<double> thermalResistancesGroutSoil2U(double chi,
                                                  double const R_ar_1,
                                                  double const R_ar_2,
                                                  double const R_g)
{
    double R_gs = compute_R_gs_2U(chi, R_g);
    double R_gg_1 = compute_R_gg_2U(chi, R_gs, R_ar_1, R_g);
    double R_gg_2 = compute_R_gg_2U(chi, R_gs, R_ar_2,
                                    R_g);  // Resulting thermal resistances.

    auto constraint = [&]() {
        return 1.0 / ((1.0 / R_gg_1) + (1.0 / (2.0 * R_gs)));
    };

    std::array<double, 3> const multiplier{chi * 0.66, chi * 0.5 * 0.66, 0.0};
    for (double m_chi : multiplier)
    {
        if (constraint() >= 0)
        {
            break;
        }
        DBUG(
            "Warning! Correction procedure was applied due to negative thermal "
            "resistance! Chi = %f.\n",
            m_chi);
        R_gs = compute_R_gs_2U(m_chi, R_g);
        R_gg_1 = compute_R_gg_2U(m_chi, R_gs, R_ar_1, R_g);
        R_gg_2 = compute_R_gg_2U(m_chi, R_gs, R_ar_2, R_g);
    }

    return {R_gg_1, R_gg_2, R_gs};
}

void BHE_2U::updateHeatTransferCoefficients(double const flow_rate)

{
    auto const tm_flow_properties = calculateThermoMechanicalFlowPropertiesPipe(
        _pipes.inlet, borehole_geometry.length, refrigerant, flow_rate);

    _flow_velocity = tm_flow_properties.velocity;
    _thermal_resistances =
        calcThermalResistances(tm_flow_properties.nusselt_number);
}

/// Nu is the Nusselt number.
std::array<double, BHE_2U::number_of_unknowns> BHE_2U::calcThermalResistances(
    double const Nu)
{
    constexpr double pi = boost::math::constants::pi<double>();

    double const lambda_r = refrigerant.thermal_conductivity;
    double const lambda_g = grout.lambda_g;
    double const lambda_p = _pipes.inlet.wall_thermal_conductivity;

    // thermal resistances due to advective flow of refrigerant in the _pipes
    // Eq. 36 in Diersch_2011_CG
    double const R_adv_i = 1.0 / (Nu * lambda_r * pi);
    double const R_adv_o = 1.0 / (Nu * lambda_r * pi);

    // thermal resistance due to thermal conductivity of the pipe wall material
    // Eq. 49
    double const R_con_a =
        std::log(_pipes.inlet.outsideDiameter() / _pipes.inlet.diameter) /
        (2.0 * pi * lambda_p);

    // the average outer diameter of the _pipes
    double const d0 = _pipes.outlet.outsideDiameter();
    double const D = borehole_geometry.diameter;
    // Eq. 38
    double const chi =
        std::log(std::sqrt(D * D + 4 * d0 * d0) / 2 / std::sqrt(2) / d0) /
        std::log(D / 2 / d0);
    // Eq. 39
    // thermal resistances of the grout
    double const R_g =
        std::acosh((D * D + d0 * d0 - 2 * _pipes.distance * _pipes.distance) /
                   (2 * D * d0)) /
        (2 * pi * lambda_g) *
        (3.098 - 4.432 * std::sqrt(2) * _pipes.distance / D +
         2.364 * 2 * _pipes.distance * _pipes.distance / D / D);
    // thermal resistance due to the grout transition.
    double const R_con_b = chi * R_g;

    // Eq. 29 and 30
    double const R_fig = 2 * R_adv_i + 2 * R_con_a + R_con_b;
    double const R_fog = 2 * R_adv_o + 2 * R_con_a + R_con_b;

    // thermal resistance due to inter-grout exchange
    double const R_ar_1 =
        std::acosh((2.0 * _pipes.distance * _pipes.distance - d0 * d0) / d0 /
                   d0) /
        (2.0 * pi * lambda_g);

    double const R_ar_2 =
        std::acosh((2.0 * 2.0 * _pipes.distance * _pipes.distance - d0 * d0) /
                   d0 / d0) /
        (2.0 * pi * lambda_g);

    std::vector<double> _intergrout_thermal_exchange;
    _intergrout_thermal_exchange =
        thermalResistancesGroutSoil2U(chi, R_ar_1, R_ar_2, R_g);

    return {{R_fig, R_fog, _intergrout_thermal_exchange[0],
             _intergrout_thermal_exchange[1], _intergrout_thermal_exchange[2]}};

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
