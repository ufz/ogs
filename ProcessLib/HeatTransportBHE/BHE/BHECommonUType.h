// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "BHECommon.h"
#include "PipeConfigurationUType.h"
#include "ThermalResistanceHelpers.h"

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
                   bool const use_python_bcs)
        : BHECommon{borehole, refrigerant, grout, flowAndTemperatureControl,
                    use_python_bcs},
          _pipes(pipes)
    {
    }

protected:
    PipeConfigurationUType const _pipes;
};
}  // namespace BHE
}  // namespace HeatTransportBHE
}  // end of namespace ProcessLib
