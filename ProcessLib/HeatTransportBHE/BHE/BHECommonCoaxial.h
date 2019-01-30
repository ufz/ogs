/**
 *
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include "BHECommon.h"
#include "PipeConfigurationCoaxial.h"

namespace ProcessLib
{
namespace HeatTransportBHE
{
namespace BHE
{
class BHECoaxialCommon : public BHECommon
{
public:
    BHECoaxialCommon(BoreholeGeometry const& borehole,
                     RefrigerantProperties const& refrigerant,
                     GroutParameters const& grout,
                     FlowAndTemperatureControl const& flowAndTemperatureControl,
                     PipeConfigurationCoaxial const& pipes)
        : BHECommon{borehole, refrigerant, grout, flowAndTemperatureControl},
          _pipes(pipes)
    {
    }
    static constexpr int number_of_unknowns = 3;
    static constexpr int number_of_grout_zones = 1;

protected:
    PipeConfigurationCoaxial const _pipes;

    /// Here we store the thermal resistances needed for computation of the heat
    /// exchange coefficients in the governing equations of BHE.
    /// These governing equations can be found in
    /// 1) Diersch (2013) FEFLOW book on page 958, M.3, or
    /// 2) Diersch (2011) Comp & Geosci 37:1122-1135, Eq. 90-97.
    std::array<double, number_of_unknowns> _thermal_resistances;

    /// Flow velocity inside the pipes and annulus. Depends on the flow_rate.
    double _flow_velocity, _flow_velocity_annulus;
};
}  // end of namespace BHE
}  // end of namespace HeatTransportBHE
}  // end of namespace ProcessLib
