/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on April 14, 2020, 3:50 PM
 */

#include "CreateRelPermNonWettingPhaseVanGenuchtenMualem.h"

#include "BaseLib/ConfigTree.h"
#include "MaterialLib/MPL/Property.h"
#include "RelPermNonWettingPhaseVanGenuchtenMualem.h"

namespace MaterialPropertyLib
{
std::unique_ptr<Property> createRelPermNonWettingPhaseVanGenuchtenMualem(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{properties__property__type}
    config.checkConfigParameter(
        "type", "RelativePermeabilityNonWettingPhaseVanGenuchtenMualem");
    DBUG("Create RelPermNonWettingPhaseVanGenuchtenMualem medium property");

    // Second access for storage.
    //! \ogs_file_param{properties__property__name}
    auto property_name = config.peekConfigParameter<std::string>("name");

    auto const residual_liquid_saturation =
        //! \ogs_file_param{properties__property__RelativePermeabilityNonWettingPhaseVanGenuchtenMualem__residual_liquid_saturation}
        config.getConfigParameter<double>("residual_liquid_saturation");
    auto const residual_gas_saturation =
        //! \ogs_file_param{properties__property__RelativePermeabilityNonWettingPhaseVanGenuchtenMualem__residual_gas_saturation}
        config.getConfigParameter<double>("residual_gas_saturation");

    auto const exponent =
        //! \ogs_file_param{properties__property__RelativePermeabilityNonWettingPhaseVanGenuchtenMualem__exponent}
        config.getConfigParameter<double>("exponent");

    auto const min_relative_permeability =
        //! \ogs_file_param{properties__property__RelativePermeabilityNonWettingPhaseVanGenuchtenMualem__min_relative_permeability}
        config.getConfigParameter<double>("min_relative_permeability");

    return std::make_unique<RelPermNonWettingPhaseVanGenuchtenMualem>(
        property_name, residual_liquid_saturation, residual_gas_saturation,
        exponent, min_relative_permeability);
}
}  // namespace MaterialPropertyLib
