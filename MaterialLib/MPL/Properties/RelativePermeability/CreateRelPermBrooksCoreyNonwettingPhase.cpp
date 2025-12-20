// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "CreateRelPermBrooksCoreyNonwettingPhase.h"

#include "BaseLib/ConfigTree.h"
#include "RelPermBrooksCoreyNonwettingPhase.h"

namespace MaterialPropertyLib
{
std::unique_ptr<RelPermBrooksCoreyNonwettingPhase>
createRelPermBrooksCoreyNonwettingPhase(BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{properties__property__type}
    config.checkConfigParameter("type", "RelPermBrooksCoreyNonwettingPhase");

    // Second access for storage.
    //! \ogs_file_param{properties__property__name}
    auto property_name = config.peekConfigParameter<std::string>("name");

    DBUG("Create RelPermBrooksCoreyNonwettingPhase medium property {:s}.",
         property_name);

    auto const residual_liquid_saturation =
        //! \ogs_file_param{properties__property__RelPermBrooksCoreyNonwettingPhase__residual_liquid_saturation}
        config.getConfigParameter<double>("residual_liquid_saturation");
    auto const residual_gas_saturation =
        //! \ogs_file_param{properties__property__RelPermBrooksCoreyNonwettingPhase__residual_gas_saturation}
        config.getConfigParameter<double>("residual_gas_saturation");
    auto const min_relative_permeability =
        //! \ogs_file_param{properties__property__RelPermBrooksCoreyNonwettingPhase__min_relative_permeability}
        config.getConfigParameter<double>("min_relative_permeability");
    auto const exponent =
        //! \ogs_file_param{properties__property__RelPermBrooksCoreyNonwettingPhase__lambda}
        config.getConfigParameter<double>("lambda");
    if (exponent <= 0.)
    {
        OGS_FATAL("Exponent 'lambda' must be positive.");
    }

    return std::make_unique<RelPermBrooksCoreyNonwettingPhase>(
        std::move(property_name),
        residual_liquid_saturation,
        residual_gas_saturation,
        min_relative_permeability,
        exponent);
}
}  // namespace MaterialPropertyLib
