// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "Pipe.h"

#include "BaseLib/ConfigTree.h"
#include "BaseLib/Error.h"

namespace ProcessLib
{
namespace HeatTransportBHE
{
namespace BHE
{
Pipe createPipe(BaseLib::ConfigTree const& config)
{
    const auto diameter =
        //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__pipes__inlet__diameter}
        config.getConfigParameter<double>("diameter");
    const auto wall_thickness =
        //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__pipes__inlet__wall_thickness}
        config.getConfigParameter<double>("wall_thickness");
    const auto wall_thermal_conductivity =
        //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__pipes__inlet__wall_thermal_conductivity}
        config.getConfigParameter<double>("wall_thermal_conductivity");

    if (diameter <= 0)
    {
        OGS_FATAL("Pipe diameter must be positive, got {:g}.", diameter);
    }
    if (wall_thickness < 0)
    {
        OGS_FATAL("Pipe wall_thickness must be non-negative, got {:g}.",
                  wall_thickness);
    }
    if (wall_thermal_conductivity <= 0)
    {
        OGS_FATAL("Pipe wall_thermal_conductivity must be positive, got {:g}.",
                  wall_thermal_conductivity);
    }

    return {diameter, wall_thickness, wall_thermal_conductivity};
}
}  // namespace BHE
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
