/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "BHE_CXC.h"
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
                 PipeConfigurationCoaxial const& pipes)
    : BHECommonCoaxial{borehole, refrigerant, grout, flowAndTemperatureControl,
                       pipes}
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

/// Nu_o is the Nusselt number of inner pipe, Nu_i is the Nusselt number of
/// annulus.
std::array<double, BHE_CXC::number_of_unknowns> BHE_CXC::calcThermalResistances(
    double const Nu_inner_pipe, double const Nu_annulus_pipe)
{
    // thermal resistances due to advective flow of refrigerant in the pipes
    auto const R_advective =
        calculateAdvectiveThermalResistance(_pipes.inner_pipe,
                                            _pipes.outer_pipe,
                                            refrigerant,
                                            Nu_inner_pipe,
                                            Nu_annulus_pipe);

    // thermal resistance due to thermal conductivity of the pipe wall material
    auto const R_conductive = calculatePipeWallThermalResistance(
        _pipes.inner_pipe, _pipes.outer_pipe);

    // thermal resistance due to the grout transition grout-soil exchange.
    auto const R = calculateGroutAndGroutSoilExchangeThermalResistance(
        _pipes.outer_pipe, grout, borehole_geometry.diameter);

    // thermal resistance due to grout-soil exchange
    double const R_gs = R.grout_soil;

    // Eq. 71 and 72
    double const R_ff = R_advective.inner_pipe_coaxial + R_advective.a_annulus +
                        R_conductive.inner_pipe_coaxial;
    double const R_fog =
        R_advective.b_annulus + R_conductive.annulus + R.conductive_b;

    return {{R_ff, R_fog, R_gs}};
}
}  // namespace BHE
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
