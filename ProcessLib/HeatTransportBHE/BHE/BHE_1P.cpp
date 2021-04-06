/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
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
    _thermal_resistances.fill(std::numeric_limits<double>::quiet_NaN());

    // Initialize thermal resistances.
    auto values = visit(
        [&](auto const& control) {
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

std::array<double, BHE_1P::number_of_unknowns> BHE_1P::pipeHeatConductions()
    const
{
    double const lambda_r = refrigerant.thermal_conductivity;
    double const rho_r = refrigerant.density;
    double const Cp_r = refrigerant.specific_heat_capacity;
    double const alpha_L = _pipe.longitudinal_dispersion_length;
    double const porosity_g = grout.porosity_g;
    double const lambda_g = grout.lambda_g;

    // Here we calculate the laplace coefficients in the governing
    // equations of the BHE.
    return {{
        // pipe, Eq. 19
        (lambda_r + rho_r * Cp_r * alpha_L * _flow_velocity),
        // grout, Eq. 21
        (1.0 - porosity_g) * lambda_g,
    }};
}

std::array<Eigen::Vector3d, BHE_1P::number_of_unknowns>
BHE_1P::pipeAdvectionVectors(Eigen::Vector3d const& elem_direction) const
{
    double const& rho_r = refrigerant.density;
    double const& Cp_r = refrigerant.specific_heat_capacity;
    Eigen::Vector3d adv_vector = rho_r * Cp_r * _flow_velocity * elem_direction;

    return {// pipe, Eq. 19
            adv_vector,
            // grout, Eq. 21
            {0, 0, 0}};
}

double BHE_1P::compute_R_gs(double const chi, double const R_g)
{
    return (1 - chi) * R_g;
}

void BHE_1P::updateHeatTransferCoefficients(double const flow_rate)

{
    auto const tm_flow_properties = calculateThermoMechanicalFlowPropertiesPipe(
        _pipe.single_pipe, borehole_geometry.length, refrigerant, flow_rate);

    _flow_velocity = tm_flow_properties.velocity;
    _thermal_resistances =
        calcThermalResistances(tm_flow_properties.nusselt_number);
}

// Nu is the Nusselt number.
std::array<double, BHE_1P::number_of_unknowns> BHE_1P::calcThermalResistances(
    double const Nu)
{
    constexpr double pi = boost::math::constants::pi<double>();

    double const lambda_r = refrigerant.thermal_conductivity;
    double const lambda_g = grout.lambda_g;
    double const lambda_p = _pipe.single_pipe.wall_thermal_conductivity;

    // thermal resistances due to advective flow of refrigerant in the pipe
    double const R_adv_i1 = 1.0 / (Nu * lambda_r * pi);

    // thermal resistance due to thermal conductivity of the pipe wall material
    double const R_con_a = std::log(_pipe.single_pipe.outsideDiameter() /
                                    _pipe.single_pipe.diameter) /
                           (2.0 * pi * lambda_p);

    // thermal resistances of the grout
    double const D = borehole_geometry.diameter;
    double const pipe_outside_diameter = _pipe.single_pipe.outsideDiameter();

    double const chi = std::log(std::sqrt(D * D + pipe_outside_diameter *
                                                      pipe_outside_diameter) /
                                std::sqrt(2) / pipe_outside_diameter) /
                       std::log(D / pipe_outside_diameter);
    double const R_g =
        std::log(D / pipe_outside_diameter) / 2 / (pi * lambda_g);

    double const R_con_b = chi * R_g;

    // thermal resistances due to grout-soil exchange
    double const R_gs = compute_R_gs(chi, R_g);

    // Eq. 29 and 30
    double const R_fg = R_adv_i1 + R_con_a + R_con_b;

    return {{R_fg, R_gs}};
}

std::array<std::pair<std::size_t /*node_id*/, int /*component*/>, 2>
BHE_1P::getBHEInflowDirichletBCNodesAndComponents(
    std::size_t const top_node_id,
    std::size_t const bottom_node_id,
    int const in_component_id)
{
    return {std::make_pair(top_node_id, in_component_id),
            std::make_pair(bottom_node_id, in_component_id)};
}

std::optional<
    std::array<std::pair<std::size_t /*node_id*/, int /*component*/>, 2>>
BHE_1P::getBHEBottomDirichletBCNodesAndComponents(
    std::size_t const /*bottom_node_id*/,
    int const /*in_component_id*/,
    int const /*out_component_id*/)
{
    return {};
}

std::array<double, BHE_1P::number_of_unknowns> BHE_1P::crossSectionAreas() const
{
    return {{_pipe.single_pipe.area(),
             borehole_geometry.area() - _pipe.single_pipe.outsideArea()}};
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
