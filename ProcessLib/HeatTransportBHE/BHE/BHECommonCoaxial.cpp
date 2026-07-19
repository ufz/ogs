// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "BHECommonCoaxial.h"

#include "BHEParameterValidation.h"
#include "ParameterLib/SpatialPosition.h"
#include "ThermalResistancesCoaxial.h"
#include "ThermoMechanicalFlowProperties.h"

namespace ProcessLib
{
namespace HeatTransportBHE
{
namespace BHE
{
BHECommonCoaxial::BHECommonCoaxial(
    BoreholeGeometry const& borehole,
    RefrigerantProperties const& refrigerant,
    GroutParameters const& grout,
    FlowAndTemperatureControl const& flowAndTemperatureControl,
    PipeConfigurationCoaxial const& pipes,
    bool const use_python_bcs)
    : BHECommon{borehole, refrigerant, grout, flowAndTemperatureControl,
                use_python_bcs},
      _pipes(pipes)
{
    cross_section_area_inner_pipe = _pipes.inner_pipe.area();
    cross_section_area_annulus =
        _pipes.outer_pipe.area() - _pipes.inner_pipe.outsideArea();
}

std::array<double, BHECommonCoaxial::number_of_unknowns>
BHECommonCoaxial::pipeHeatCapacities() const
{
    double const rho_r = refrigerant.density;
    double const specific_heat_capacity = refrigerant.specific_heat_capacity;
    double const rho_g = grout.rho_g;
    double const porosity_g = grout.porosity_g;
    double const heat_cap_g = grout.heat_cap_g;

    return {{/*i*/ rho_r * specific_heat_capacity,
             /*o*/ rho_r * specific_heat_capacity,
             /*g*/ (1.0 - porosity_g) * rho_g * heat_cap_g}};
}

double BHECommonCoaxial::updateFlowRateAndTemperature(double const T_out,
                                                      double const current_time)
{
    auto values =
        visit([&](auto const& control) { return control(T_out, current_time); },
              flowAndTemperatureControl);
    updateHeatTransferCoefficients(values.flow_rate);
    return values.temperature;
}

std::array<double, BHECommonCoaxial::number_of_unknowns>
BHECommonCoaxial::pipeHeatConductions(
    ParameterLib::SpatialPosition const& /*pos*/) const
{
    double const lambda_r = refrigerant.thermal_conductivity;
    double const rho_r = refrigerant.density;
    double const Cp_r = refrigerant.specific_heat_capacity;
    double const alpha_L = _pipes.longitudinal_dispersion_length;
    double const porosity_g = grout.porosity_g;
    double const lambda_g = grout.lambda_g;

    auto const legs = flowLegs();
    return {{(lambda_r + rho_r * Cp_r * alpha_L * std::abs(legs[0])),
             (lambda_r + rho_r * Cp_r * alpha_L * std::abs(legs[1])),
             (1.0 - porosity_g) * lambda_g}};
}

std::array<Eigen::Vector3d, BHECommonCoaxial::number_of_unknowns>
BHECommonCoaxial::pipeAdvectionVectors(
    Eigen::Vector3d const& elem_direction) const
{
    double const rho_r = refrigerant.density;
    double const Cp_r = refrigerant.specific_heat_capacity;

    auto const legs = flowLegs();
    auto leg_adv = [&](double const v_signed) -> Eigen::Vector3d
    { return rho_r * Cp_r * v_signed * elem_direction; };

    return {leg_adv(legs[0]), leg_adv(legs[1]), {0, 0, 0}};
}

std::vector<double> BHECommonCoaxial::calcThermalResistances(
    double const Nu_inner_pipe, double const Nu_annulus_pipe,
    ParameterLib::SpatialPosition const& pos) const
{
    // thermal resistances due to advective flow of refrigerant in the pipes
    auto const R_advective = calculateAdvectiveThermalResistance(
        _pipes.inner_pipe, _pipes.outer_pipe, refrigerant, Nu_inner_pipe,
        Nu_annulus_pipe);

    // thermal resistance due to thermal conductivity of the pipe wall material
    // t=0.0: borehole properties are physically time-invariant; genuinely
    // time-varying parameter types are rejected in createPipe.
    double const lambda_p_inner =
        sampleStrictPositive(_pipes.inner_pipe.wall_thermal_conductivity, 0.0,
                             pos, "wall_thermal_conductivity (inner pipe)");
    double const lambda_p_outer =
        sampleStrictPositive(_pipes.outer_pipe.wall_thermal_conductivity, 0.0,
                             pos, "wall_thermal_conductivity (outer pipe)");
    auto const R_conductive = calculatePipeWallThermalResistance(
        _pipes.inner_pipe, lambda_p_inner, _pipes.outer_pipe, lambda_p_outer);

    // thermal resistance due to the grout transition and grout-soil exchange.
    double const D = sampleStrictPositive(borehole_geometry.diameter, 0.0, pos,
                                          "borehole_diameter");
    checkBoreholeVsPipeDiameter(D, _pipes.outer_pipe.outsideDiameter(), pos,
                                "coaxial grout resistance");
    auto const R = calculateGroutAndGroutSoilExchangeThermalResistance(
        _pipes.outer_pipe, grout, D);

    double const R_gs = R.grout_soil;

    double const R_ff = R_advective.inner_pipe_coaxial + R_advective.a_annulus +
                        R_conductive.inner_pipe_coaxial;
    double const R_fg =
        R_advective.b_annulus + R_conductive.annulus + R.conductive_b;

    return getThermalResistances(R_gs, R_ff, R_fg);
}

std::vector<double> BHECommonCoaxial::thermalResistances(
    ParameterLib::SpatialPosition const& pos) const
{
    return calcThermalResistances(cached_nu_inner, cached_nu_annulus, pos);
}

std::array<std::pair<std::size_t, int>, 2>
BHECommonCoaxial::getBHEInflowDirichletBCNodesAndComponents(
    std::size_t const top_node_id,
    std::size_t const /*bottom_node_id*/,
    int const in_component_id)
{
    return {std::make_pair(top_node_id, in_component_id),
            std::make_pair(top_node_id, in_component_id + 1)};
}

std::optional<std::array<std::pair<std::size_t, int>, 2>>
BHECommonCoaxial::getBHEBottomDirichletBCNodesAndComponents(
    std::size_t const bottom_node_id, int const in_component_id,
    int const out_component_id)
{
    return {{std::make_pair(bottom_node_id, in_component_id),
             std::make_pair(bottom_node_id, out_component_id)}};
}

void BHECommonCoaxial::updateHeatTransferCoefficients(double const flow_rate)
{
    auto const tm_flow_inner = calculateThermoMechanicalFlowPropertiesPipe(
        _pipes.inner_pipe, borehole_geometry.length, refrigerant, flow_rate);

    auto const tm_flow_annulus = calculateThermoMechanicalFlowPropertiesAnnulus(
        _pipes.inner_pipe, _pipes.outer_pipe, borehole_geometry.length,
        refrigerant, flow_rate);

    cached_nu_inner = tm_flow_inner.nusselt_number;
    cached_nu_annulus = tm_flow_annulus.nusselt_number;
    assignVelocities(tm_flow_inner.velocity, tm_flow_annulus.velocity);
}
}  // namespace BHE
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
