/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "BHE_CXC.h"

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
BHE_CXC::BHE_CXC(BoreholeGeometry const& borehole,
                 RefrigerantProperties const& refrigerant,
                 GroutParameters const& grout,
                 FlowAndTemperatureControl const& flowAndTemperatureControl,
                 PipeConfigurationCXC const& pipes)
    : BHECommon{borehole, refrigerant, grout, flowAndTemperatureControl},
      _pipes(pipes)
{
    // Initialize thermal resistances.
    auto values = apply_visitor(
        [&](auto const& control) {
            return control(refrigerant.reference_temperature,
                           0. /* initial time */);
        },
        flowAndTemperatureControl);
    updateHeatTransferCoefficients(values.flow_rate);
}

std::array<double, BHE_CXC::number_of_unknowns> BHE_CXC::pipeHeatCapacities()
    const
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

std::array<double, BHE_CXC::number_of_unknowns> BHE_CXC::pipeHeatConductions()
    const
{
    double const& lambda_r = refrigerant.thermal_conductivity;
    double const& rho_r = refrigerant.density;
    double const& Cp_r = refrigerant.specific_heat_capacity;
    double const& alpha_L = _pipes.longitudinal_dispersion_length;
    double const& porosity_g = grout.porosity_g;
    double const& lambda_g = grout.lambda_g;

    double const velocity_norm_inner_pipe = std::abs(_flow_velocity);
    double const velocity_norm_annulus = std::abs(_flow_velocity_annulus);

    // Here we calculate the laplace coefficients in the governing
    // equations of BHE. These governing equations can be found in
    // 1) Diersch (2013) FEFLOW book on page 952, M.120-122, or
    // 2) Diersch (2011) Comp & Geosci 37:1122-1135, Eq. 26-28.
    return {{// pipe i, Eq. 26
             (lambda_r + rho_r * Cp_r * alpha_L * velocity_norm_inner_pipe),
             // pipe o, Eq. 27
             (lambda_r + rho_r * Cp_r * alpha_L * velocity_norm_annulus),
             // pipe g, Eq. 28
             (1.0 - porosity_g) * lambda_g}};
}

std::array<Eigen::Vector3d, BHE_CXC::number_of_unknowns>
BHE_CXC::pipeAdvectionVectors() const
{
    double const& rho_r = refrigerant.density;
    double const& Cp_r = refrigerant.specific_heat_capacity;

    return {{// pipe i, Eq. 26
             {0, 0, -rho_r * Cp_r * _flow_velocity},
             // pipe o, Eq. 27
             {0, 0, rho_r * Cp_r * _flow_velocity_annulus},
             // grout g, Eq. 28
             {0, 0, 0}}};
}

constexpr std::pair<int, int> BHE_CXC::inflow_outflow_bc_component_ids[];

void BHE_CXC::updateHeatTransferCoefficients(double const flow_rate)

{
    auto const tm_flow_properties_annulus =
        calculateThermoMechanicalFlowPropertiesAnnulus(_pipes.inner_inflow_pipe,
                                                       _pipes.outer_pipe,
                                                       borehole_geometry.length,
                                                       refrigerant,
                                                       flow_rate);

    _flow_velocity_annulus = tm_flow_properties_annulus.velocity;

    auto const tm_flow_properties =
        calculateThermoMechanicalFlowPropertiesPipe(_pipes.inner_inflow_pipe,
                                                    borehole_geometry.length,
                                                    refrigerant,
                                                    flow_rate);

    _flow_velocity = tm_flow_properties.velocity;
    _thermal_resistances =
        calcThermalResistances(tm_flow_properties_annulus.nusselt_number,
                               tm_flow_properties.nusselt_number);
}

/// Nu_o is the Nusselt number of inner pipe, Nu_i is the Nusselt number of
/// annulus.
std::array<double, BHE_CXC::number_of_unknowns> BHE_CXC::calcThermalResistances(
    double const Nu_annulus, double const Nu_inner_inflow_pipe)
{
    static constexpr double pi = boost::math::constants::pi<double>();

    double const& D = borehole_geometry.diameter;
    double const& lambda_r = refrigerant.thermal_conductivity;
    double const& lambda_g = grout.lambda_g;
    double const& d_o = _pipes.outer_pipe.diameter;
    double const& d_i = _pipes.inner_inflow_pipe.diameter;
    double const& b_o = _pipes.outer_pipe.wall_thickness;
    double const& b_i = _pipes.inner_inflow_pipe.wall_thickness;
    double const& lambda_p_o = _pipes.outer_pipe.wall_thermal_conductivity;
    double const& lambda_p_i =
        _pipes.inner_inflow_pipe.wall_thermal_conductivity;

    // thermal resistances due to advective flow of refrigerant in the _pipes
    // Eq. 73, 74, 75 in Diersch_2011_CG
    double const d_h = d_o - d_i - 2 * b_i;
    double const d_i_o = d_i + 2 * b_i;
    double const d_o_o = d_o + 2 * b_o;
    double const R_adv_i = 1.0 / (Nu_inner_inflow_pipe * lambda_r * pi);
    double const R_adv_a_o = 1.0 / (Nu_annulus * lambda_r * pi) * (d_h / d_i_o);
    double const R_adv_b_o = 1.0 / (Nu_annulus * lambda_r * pi) * (d_h / d_o);

    // thermal resistance due to thermal conductivity of the pipe wall material
    // Eq. 81
    double const R_con_o = std::log(d_o_o / d_o) / (2.0 * pi * lambda_p_o);
    double const R_con_i = std::log(d_i_o / d_i) / (2.0 * pi * lambda_p_i);

    // Eq. 83
    double const chi =
        std::log(std::sqrt(D * D + d_o_o * d_o_o) / std::sqrt(2) / d_o_o) /
        std::log(D / d_o_o);
    // Eq. 84
    // thermal resistances of the grout
    double const R_g = std::log(D / d_o_o) / 2 / pi / lambda_g;
    // Eq. 82
    // thermal resistance due to the grout transition.
    double const R_con_b = chi * R_g;
    // Eq. 71 and 72
    double const R_ff = R_adv_i + R_adv_a_o + R_con_i;
    double const R_fog = R_adv_b_o + R_con_o + R_con_b;

    // thermal resistance due to grout-soil exchange
    double const R_gs = (1 - chi) * R_g;

    return {{R_ff, R_fog, R_gs}};

    // keep the following lines------------------------------------------------
    // when debuffing the code, printing the R and phi values are needed--------
    // std::cout << "Rfog =" << R_fog << " Rff =" << R_ff << " Rgs =" << R_gs
    // << "\n"; double phi_fog = 1.0 / (R_fog * S_i);
    // double phi_ff = 1.0 / (R_ff * S_g1); double phi_gs = 1.0 / (R_gs * S_gs);
    // std::cout << "phi_fog ="
    // << phi_fog << " phi_ff =" << phi_ff << " phi_gs =" << phi_gs << "\n";
    // -------------------------------------------------------------------------
}

double BHE_CXC::updateFlowRateAndTemperature(double const T_out,
                                             double const current_time)
{
    auto values = apply_visitor(
        [&](auto const& control) { return control(T_out, current_time); },
        flowAndTemperatureControl);
    updateHeatTransferCoefficients(values.flow_rate);
    return values.temperature;
}
}  // namespace BHE
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
