/**
 * \file
 * \author Norbert Grunwald
 * \date   01.12.2020
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "BaseLib/ConfigTree.h"
#include "RelPermUdell.h"

namespace MaterialPropertyLib
{
std::unique_ptr<RelPermUdell> createRelPermUdell(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{properties__property__type}
    config.checkConfigParameter("type", "RelativePermeabilityUdell");

    // Second access for storage.
    //! \ogs_file_param{properties__property__name}
    auto property_name = config.peekConfigParameter<std::string>("name");

    DBUG("Create RelPermUdell medium property {:s}.", property_name);

    auto const residual_liquid_saturation =
        //! \ogs_file_param{properties__property__RelPermUdell__residual_liquid_saturation}
        config.getConfigParameter<double>("residual_liquid_saturation");
    auto const residual_gas_saturation =
        //! \ogs_file_param{properties__property__RelPermUdell__residual_gas_saturation}
        config.getConfigParameter<double>("residual_gas_saturation");
    auto const min_relative_permeability_liquid =
        //! \ogs_file_param{properties__property__RelPermUdell__min_relative_permeability_liquid}
        config.getConfigParameter<double>("min_relative_permeability_liquid");
    auto const min_relative_permeability_gas =
        //! \ogs_file_param{properties__property__RelPermUdell__min_relative_permeability_gas}
        config.getConfigParameter<double>("min_relative_permeability_gas");

    if ((min_relative_permeability_liquid < 0) ||
        (min_relative_permeability_gas < 0))
    {
        OGS_FATAL("Minimal relative permeabilities must be non-negative.");
    }

    return std::make_unique<RelPermUdell>(
        std::move(property_name), residual_liquid_saturation,
        residual_gas_saturation, min_relative_permeability_liquid,
        min_relative_permeability_gas);
}
}  // namespace MaterialPropertyLib
