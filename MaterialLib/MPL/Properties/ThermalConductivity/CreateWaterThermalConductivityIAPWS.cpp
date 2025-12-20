// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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
