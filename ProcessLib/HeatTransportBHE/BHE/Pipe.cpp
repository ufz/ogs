// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "Pipe.h"

#include "BHEParameterValidation.h"
#include "BaseLib/ConfigTree.h"
#include "BaseLib/Error.h"
#include "ParameterLib/Parameter.h"

namespace ProcessLib
{
namespace HeatTransportBHE
{
namespace BHE
{
Pipe createPipe(
    BaseLib::ConfigTree const& config,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>>& parameters)
{
    const auto diameter =
        //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__pipes__inlet__diameter}
        config.getConfigParameter<double>("diameter");
    const auto wall_thickness =
        //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__pipes__inlet__wall_thickness}
        config.getConfigParameter<double>("wall_thickness");

    auto const wall_thermal_conductivity_str =
        //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__pipes__inlet__wall_thermal_conductivity}
        config.getConfigParameter<std::string>("wall_thermal_conductivity");
    auto const& wall_thermal_conductivity_param =
        ParameterLib::getNamedOrCreateInlineParameter(
            wall_thermal_conductivity_str, parameters,
            "pipe_wall_thermal_conductivity", "pipe_wall_thermal_conductivity");

    if (diameter <= 0)
    {
        OGS_FATAL("Pipe diameter must be positive, got {:g}.", diameter);
    }
    if (wall_thickness < 0)
    {
        OGS_FATAL("Pipe wall_thickness must be non-negative, got {:g}.",
                  wall_thickness);
    }

    // Borehole pipe wall conductivity is physically time-invariant (fixed at
    // drilling / installation).  Reject time-varying / mesh-node / non-scalar
    // parameter types; warn for FunctionParameter (may be spatial-only).
    validateScalarTimeInvariantBHEParameter(wall_thermal_conductivity_param,
                                            "wall_thermal_conductivity");

    return {diameter, wall_thickness, wall_thermal_conductivity_param};
}
}  // namespace BHE
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
