/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "BaseLib/ConfigTree.h"
#include "RelPermVanGenuchten.h"

namespace MaterialPropertyLib
{
std::unique_ptr<RelPermVanGenuchten> createRelPermVanGenuchten(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{properties__property__type}
    config.checkConfigParameter("type", "RelativePermeabilityVanGenuchten");

    // Second access for storage.
    //! \ogs_file_param{properties__property__name}
    auto property_name = config.peekConfigParameter<std::string>("name");

    DBUG("Create RelativePermeabilityVanGenuchten medium property {:s}.",
         property_name);

    auto const residual_liquid_saturation =
        //! \ogs_file_param{properties__property__RelativePermeabilityVanGenuchten__residual_liquid_saturation}
        config.getConfigParameter<double>("residual_liquid_saturation");
    auto const residual_gas_saturation =
        //! \ogs_file_param{properties__property__RelativePermeabilityVanGenuchten__residual_gas_saturation}
        config.getConfigParameter<double>("residual_gas_saturation");
    auto const min_relative_permeability_liquid =
        //! \ogs_file_param{properties__property__RelativePermeabilityVanGenuchten__minimum_relative_permeability_liquid}
        config.getConfigParameter<double>(
            "minimum_relative_permeability_liquid");
    auto const exponent =
        //! \ogs_file_param{properties__property__RelativePermeabilityVanGenuchten__exponent}
        config.getConfigParameter<double>("exponent");
    if (exponent <= 0. || exponent >= 1.)
    {
        OGS_FATAL("Exponent must be in the (0, 1) range.");
    }

    return std::make_unique<RelPermVanGenuchten>(
        std::move(property_name),
        residual_liquid_saturation,
        residual_gas_saturation,
        min_relative_permeability_liquid,
        exponent);
}
}  // namespace MaterialPropertyLib
