// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "BoreholeGeometry.h"

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
BoreholeGeometry createBoreholeGeometry(
    BaseLib::ConfigTree const& config,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>>& parameters)
{
    const auto borehole_length =
        //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__borehole__length}
        config.getConfigParameter<double>("length");
    if (borehole_length <= 0)
    {
        OGS_FATAL("Borehole length must be positive, got {:g}.",
                  borehole_length);
    }

    auto const diameter_parameter_or_value =
        //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__borehole__diameter}
        config.getConfigParameter<std::string>("diameter");
    auto const& diameter_parameter =
        ParameterLib::getNamedOrCreateInlineParameter(
            diameter_parameter_or_value, parameters, "borehole_diameter",
            "borehole_diameter");

    // Borehole diameter is physically time-invariant (fixed at drilling).
    // Reject time-varying / mesh-node / non-scalar parameter types; warn for
    // FunctionParameter (which may be spatial-only).
    validateScalarTimeInvariantBHEParameter(diameter_parameter,
                                            "borehole_diameter");

    return {borehole_length, diameter_parameter};
}
}  // namespace BHE
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
