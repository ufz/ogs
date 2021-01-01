/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on March 20, 2020, 1:31 PM
 */

#include "CreateCapillaryPressureVanGenuchten.h"

#include "BaseLib/ConfigTree.h"
#include "CapillaryPressureVanGenuchten.h"
#include "MaterialLib/MPL/Property.h"

namespace MaterialPropertyLib
{
std::unique_ptr<Property> createCapillaryPressureVanGenuchten(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{properties__property__type}
    config.checkConfigParameter("type", "CapillaryPressureVanGenuchten");

    // Second access for storage.
    //! \ogs_file_param{properties__property__name}
    auto property_name = config.peekConfigParameter<std::string>("name");

    DBUG("Create CapillaryPressureVanGenuchten medium property {:s}.",
         property_name);

    auto const residual_liquid_saturation =
        //! \ogs_file_param{properties__property__CapillaryPressureVanGenuchten__residual_liquid_saturation}
        config.getConfigParameter<double>("residual_liquid_saturation");
    auto const residual_gas_saturation =
        //! \ogs_file_param{properties__property__CapillaryPressureVanGenuchten__residual_gas_saturation}
        config.getConfigParameter<double>("residual_gas_saturation");
    auto const exponent =
        //! \ogs_file_param{properties__property__CapillaryPressureVanGenuchten__exponent}
        config.getConfigParameter<double>("exponent");
    auto const p_b =
        //! \ogs_file_param{properties__property__CapillaryPressureVanGenuchten__p_b}
        config.getConfigParameter<double>("p_b");
    auto const maximum_capillary_pressure =
        //! \ogs_file_param{properties__property__CapillaryPressureVanGenuchten__maximum_capillary_pressure}
        config.getConfigParameter<double>("maximum_capillary_pressure");

    return std::make_unique<CapillaryPressureVanGenuchten>(
        std::move(property_name), residual_liquid_saturation,
        residual_gas_saturation, exponent, p_b, maximum_capillary_pressure);
}
}  // namespace MaterialPropertyLib
