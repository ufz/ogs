/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
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
std::unique_ptr<Property>
    createCapillaryPressureVanGenuchten(BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{properties__property__type}
    config.checkConfigParameter("type", "CapillaryPressureVanGenuchten");

    DBUG("Create CapillaryPressureVanGenuchten medium property");

    auto const residual_liquid_saturation =
        //! \ogs_file_param{properties__property__CapillaryPressureVanGenuchten__residual_liquid_saturation}
        config.getConfigParameter<double>("residual_liquid_saturation");
    auto const residual_gas_saturation =
        //! \ogs_file_param{properties__property__CapillaryPressureVanGenuchten__residual_gas_saturation}
        config.getConfigParameter<double>("residual_gas_saturation");
    auto const exponent =
        //! \ogs_file_param{properties__property__CapillaryPressureVanGenuchten__exponent}
        config.getConfigParameter<double>("exponent");
    auto const entry_pressure =
        //! \ogs_file_param{properties__property__CapillaryPressureVanGenuchten__entry_pressure}
        config.getConfigParameter<double>("entry_pressure");
    auto const maximum_capillary_pressure =
        //! \ogs_file_param{properties__property__CapillaryPressureVanGenuchten__maximum_capillary_pressure}
        config.getConfigParameter<double>("maximum_capillary_pressure");

    auto const has_regularization =
        //! \ogs_file_param{properties__property__CapillaryPressureVanGenuchten__regularization}
        config.getConfigParameter<bool>("regularization", false);

    if (has_regularization)
    {
        return std::make_unique<CapillaryPressureRegularizedVanGenuchten>(
            residual_liquid_saturation, residual_gas_saturation, exponent,
            entry_pressure, maximum_capillary_pressure);
    }

    return std::make_unique<CapillaryPressureVanGenuchten>(
        residual_liquid_saturation, residual_gas_saturation, exponent,
        entry_pressure, maximum_capillary_pressure);
}
}  // namespace MaterialPropertyLib
