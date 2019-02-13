/**
 *
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <Eigen/Eigen>
#include "BHECommon.h"
#include "FlowAndTemperatureControl.h"
#include "Physics.h"
#include "PipeConfigurationCoaxial.h"
#include "ThermalResistancesCoaxial.h"
#include "ThermoMechanicalFlowProperties.h"

namespace ProcessLib
{
namespace HeatTransportBHE
{
namespace BHE
{
class BHECommonCoaxial : public BHECommon
{
public:
    BHECommonCoaxial(BoreholeGeometry const& borehole,
                     RefrigerantProperties const& refrigerant,
                     GroutParameters const& grout,
                     FlowAndTemperatureControl const& flowAndTemperatureControl,
                     PipeConfigurationCoaxial const& pipes)
        : BHECommon{borehole, refrigerant, grout, flowAndTemperatureControl},
          _pipes(pipes)
    {
        cross_section_area_inner_pipe = _pipes.inner_pipe.area();
        cross_section_area_annulus =
            _pipes.outer_pipe.area() - _pipes.inner_pipe.outsideArea();
        cross_section_area_grout =
            borehole_geometry.area() - _pipes.outer_pipe.outsideArea();
    }

    static constexpr int number_of_unknowns = 3;
    static constexpr int number_of_grout_zones = 1;

    double thermalResistance(int const unknown_index) const
    {
        return _thermal_resistances[unknown_index];
    }

    double updateFlowRateAndTemperature(double T_out, double current_time);

    std::array<double, number_of_unknowns> calcThermalResistances(
        double const Nu_inner_pipe, double const Nu_annulus_pipe);

    std::array<double, number_of_unknowns> pipeHeatCapacities() const;

    static constexpr std::pair<int, int> inflow_outflow_bc_component_ids[] = {
        {0, 1}};

    std::array<double, number_of_unknowns> pipeHeatConductions() const;

    std::array<Eigen::Vector3d, number_of_unknowns> pipeAdvectionVectors()
        const;

    double cross_section_area_inner_pipe, cross_section_area_annulus,
        cross_section_area_grout;

protected:
    void updateHeatTransferCoefficients(double const flow_rate)
    {
        auto const tm_flow_properties_annulus =
            calculateThermoMechanicalFlowPropertiesAnnulus(
                _pipes.inner_pipe,
                _pipes.outer_pipe,
                borehole_geometry.length,
                refrigerant,
                flow_rate);

        _flow_velocity_annulus = tm_flow_properties_annulus.velocity;

        auto const tm_flow_properties =
            calculateThermoMechanicalFlowPropertiesPipe(
                _pipes.inner_pipe,
                borehole_geometry.length,
                refrigerant,
                flow_rate);

        _flow_velocity_inner = tm_flow_properties.velocity;

        _thermal_resistances =
            calcThermalResistances(tm_flow_properties.nusselt_number,
                                   tm_flow_properties_annulus.nusselt_number);
    }

    PipeConfigurationCoaxial const _pipes;

    virtual std::array<double, 2> velocities() const = 0;

    virtual std::array<double, number_of_unknowns> getThermalResistances(
        double const& R_gs, double const& R_ff, double const& R_fg) const = 0;

    /// Here we store the thermal resistances needed for computation of the heat
    /// exchange coefficients in the governing equations of BHE.
    /// These governing equations can be found in
    /// 1) Diersch (2013) FEFLOW book on page 958, M.3, or
    /// 2) Diersch (2011) Comp & Geosci 37:1122-1135, Eq. 90-97.
    std::array<double, number_of_unknowns> _thermal_resistances;

    /// Flow velocity inside the pipes and annulus. Depends on the flow_rate.
    double _flow_velocity_inner, _flow_velocity_annulus;
};
}  // namespace BHE
}  // namespace HeatTransportBHE
}  // end of namespace ProcessLib
