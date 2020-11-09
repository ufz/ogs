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
#include "SaturationExponential.h"

namespace MaterialPropertyLib
{
std::unique_ptr<SaturationExponential> createSaturationExponential(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{properties__property__type}
    config.checkConfigParameter("type", "SaturationExponential");

    // Second access for storage.
    //! \ogs_file_param{properties__property__name}
    auto property_name = config.peekConfigParameter<std::string>("name");

    DBUG("Create SaturationExponential medium property {:s}.", property_name);

    auto const residual_liquid_saturation =
        //! \ogs_file_param{properties__property__SaturationExponential__residual_liquid_saturation}
        config.getConfigParameter<double>("residual_liquid_saturation");
    auto const residual_gas_saturation =
        //! \ogs_file_param{properties__property__SaturationExponential__residual_gas_saturation}
        config.getConfigParameter<double>("residual_gas_saturation");
    auto const p_cap_ref =
        //! \ogs_file_param{properties__property__SaturationExponential__reference_capillary_pressure}
        config.getConfigParameter<double>("reference_capillary_pressure");
    auto const exponent =
        //! \ogs_file_param{properties__property__SaturationExponential__exponent}
        config.getConfigParameter<double>("exponent");

    return std::make_unique<SaturationExponential>(
        std::move(property_name), residual_liquid_saturation,
        residual_gas_saturation, p_cap_ref, exponent);
}
}  // namespace MaterialPropertyLib
