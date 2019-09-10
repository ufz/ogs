/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <Eigen/Eigen>
#include "BHECommon.h"
#include "FlowAndTemperatureControl.h"
#include "PipeConfigurationUType.h"

namespace ProcessLib
{
namespace HeatTransportBHE
{
namespace BHE
{
class BHECommonUType : public BHECommon
{
public:
    BHECommonUType(BoreholeGeometry const& borehole,
                   RefrigerantProperties const& refrigerant,
                   GroutParameters const& grout,
                   FlowAndTemperatureControl const& flowAndTemperatureControl,
                   PipeConfigurationUType const& pipes,
                   bool const ifUsePythonBC)
        : BHECommon{borehole, refrigerant, grout, flowAndTemperatureControl,
                    ifUsePythonBC},
          _pipes(pipes)
    {
    }

protected:
    PipeConfigurationUType const _pipes;

    /// Flow velocity inside the pipes. Depends on the flow_rate.
    double _flow_velocity = std::numeric_limits<double>::quiet_NaN();
};
}  // namespace BHE
}  // namespace HeatTransportBHE
}  // end of namespace ProcessLib
