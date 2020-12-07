/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "BaseLib/ConfigTree.h"
#include "SaturationDependentHeatConduction.h"

namespace MaterialPropertyLib
{
std::unique_ptr<SaturationDependentHeatConduction>
createSaturationDependentHeatConduction(BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{properties__property__type}
    config.checkConfigParameter("type", "SaturationDependentHeatConduction");

    // Second access for storage.
    //! \ogs_file_param{properties__property__name}
    auto property_name = config.peekConfigParameter<std::string>("name");

    DBUG("Create saturation dependent thermal_conductivity property {:s}.",
         property_name);

    auto const K_dry =
        //! \ogs_file_param{properties__property__SaturationDependentHeatConduction__dry}
        config.getConfigParameter<double>("dry");

    auto const K_wet =
        //! \ogs_file_param{properties__property__SaturationDependentHeatConduction__wet}
        config.getConfigParameter<double>("wet");

    return std::make_unique<
        MaterialPropertyLib::SaturationDependentHeatConduction>(
        std::move(property_name), K_dry, K_wet);
}
}  // namespace MaterialPropertyLib
