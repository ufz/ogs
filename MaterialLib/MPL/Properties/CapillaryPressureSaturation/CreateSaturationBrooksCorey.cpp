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
#include "SaturationBrooksCorey.h"

namespace MaterialPropertyLib
{
std::unique_ptr<SaturationBrooksCorey> createSaturationBrooksCorey(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{properties__property__type}
    config.checkConfigParameter("type", "SaturationBrooksCorey");

    // Second access for storage.
    //! \ogs_file_param{properties__property__name}
    auto property_name = config.peekConfigParameter<std::string>("name");

    DBUG("Create SaturationBrooksCorey medium property {:s}.", property_name);

    auto const residual_liquid_saturation =
        //! \ogs_file_param{properties__property__SaturationBrooksCorey__residual_liquid_saturation}
        config.getConfigParameter<double>("residual_liquid_saturation");
    auto const residual_gas_saturation =
        //! \ogs_file_param{properties__property__SaturationBrooksCorey__residual_gas_saturation}
        config.getConfigParameter<double>("residual_gas_saturation");
    auto const exponent =
        //! \ogs_file_param{properties__property__SaturationBrooksCorey__lambda}
        config.getConfigParameter<double>("lambda");
    auto const entry_pressure =
        //! \ogs_file_param{properties__property__SaturationBrooksCorey__entry_pressure}
        config.getConfigParameter<double>("entry_pressure");

    return std::make_unique<SaturationBrooksCorey>(
        std::move(property_name), residual_liquid_saturation,
        residual_gas_saturation, exponent, entry_pressure);
}
}  // namespace MaterialPropertyLib
