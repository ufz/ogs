// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "CreateWaterVapourLatentHeatWithCriticalTemperature.h"

#include "BaseLib/ConfigTree.h"
#include "MaterialLib/MPL/Property.h"
#include "WaterVapourLatentHeatWithCriticalTemperature.h"

namespace MaterialPropertyLib
{
std::unique_ptr<Property> createWaterVapourLatentHeatWithCriticalTemperature(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{properties__property__type}
    config.checkConfigParameter("type",
                                "WaterVapourLatentHeatWithCriticalTemperature");
    DBUG("Create WaterVapourLatentHeatWithCriticalTemperature phase property");

    // Second access for storage.
    //! \ogs_file_param{properties__property__name}
    auto property_name = config.peekConfigParameter<std::string>("name");

    //! \ogs_file_param_special{properties__property__WaterVapourLatentHeatWithCriticalTemperature}
    return std::make_unique<WaterVapourLatentHeatWithCriticalTemperature>(
        std::move(property_name));
}
}  // namespace MaterialPropertyLib
