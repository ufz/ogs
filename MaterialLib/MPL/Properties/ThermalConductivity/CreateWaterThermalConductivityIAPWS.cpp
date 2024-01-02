/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on March 4, 2021, 4:38 PM
 */

#include "CreateWaterThermalConductivityIAPWS.h"

#include "BaseLib/ConfigTree.h"
#include "MaterialLib/MPL/Property.h"
#include "WaterThermalConductivityIAPWS.h"

namespace MaterialPropertyLib
{
std::unique_ptr<Property> createWaterThermalConductivityIAPWS(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{properties__property__type}
    config.checkConfigParameter("type", "WaterThermalConductivityIAPWS");
    DBUG("Create WaterThermalConductivityIAPWS phase property");

    // Second access for storage.
    //! \ogs_file_param{properties__property__name}
    auto property_name = config.peekConfigParameter<std::string>("name");

    //! \ogs_file_param_special{properties__property__WaterThermalConductivityIAPWS}
    return std::make_unique<WaterThermalConductivityIAPWS>(
        std::move(property_name));
}
}  // namespace MaterialPropertyLib
