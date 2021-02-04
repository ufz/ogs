/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on January 26, 2021, 11:11 AM
 */

#include "BaseLib/ConfigTree.h"
#include "RelPermUdellNonwettingPhase.h"

namespace MaterialPropertyLib
{
std::unique_ptr<RelPermUdellNonwettingPhase> createRelPermUdellNonwettingPhase(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{properties__property__type}
    config.checkConfigParameter("type",
                                "RelativePermeabilityUdellNonwettingPhase");

    // Second access for storage.
    //! \ogs_file_param{properties__property__name}
    auto property_name = config.peekConfigParameter<std::string>("name");

    DBUG("Create RelPermUdellNonwettingPhase medium property {:s}.",
         property_name);

    auto const residual_liquid_saturation =
        //! \ogs_file_param{properties__property__RelativePermeabilityUdellNonwettingPhase__residual_liquid_saturation}
        config.getConfigParameter<double>("residual_liquid_saturation");
    auto const residual_gas_saturation =
        //! \ogs_file_param{properties__property__RelativePermeabilityUdellNonwettingPhase__residual_gas_saturation}
        config.getConfigParameter<double>("residual_gas_saturation");
    auto const min_relative_permeability =
        //! \ogs_file_param{properties__property__RelativePermeabilityUdellNonwettingPhase__min_relative_permeability}
        config.getConfigParameter<double>("min_relative_permeability");

    if (min_relative_permeability < 0)
    {
        OGS_FATAL("Minimal relative permeability must be non-negative.");
    }

    return std::make_unique<RelPermUdellNonwettingPhase>(
        std::move(property_name), residual_liquid_saturation,
        residual_gas_saturation, min_relative_permeability);
}
}  // namespace MaterialPropertyLib
