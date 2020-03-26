/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "BaseLib/ConfigTree.h"
#include "SaturationVanGenuchten.h"

namespace MaterialPropertyLib
{
std::unique_ptr<SaturationVanGenuchten> createSaturationVanGenuchten(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{properties__property__type}
    config.checkConfigParameter("type", "SaturationVanGenuchten");

    DBUG("Create SaturationVanGenuchten medium property");

    auto const residual_liquid_saturation =
        //! \ogs_file_param{properties__property__SaturationVanGenuchten__residual_liquid_saturation}
        config.getConfigParameter<double>("residual_liquid_saturation");
    auto const residual_gas_saturation =
        //! \ogs_file_param{properties__property__SaturationVanGenuchten__residual_gas_saturation}
        config.getConfigParameter<double>("residual_gas_saturation");
    auto const exponent =
        //! \ogs_file_param{properties__property__SaturationVanGenuchten__exponent}
        config.getConfigParameter<double>("exponent");
    auto const entry_pressure =
        //! \ogs_file_param{properties__property__SaturationVanGenuchten__entry_pressure}
        config.getConfigParameter<double>("entry_pressure");

    return std::make_unique<SaturationVanGenuchten>(residual_liquid_saturation,
                                                    residual_gas_saturation,
                                                    exponent, entry_pressure);
}
}  // namespace MaterialPropertyLib
