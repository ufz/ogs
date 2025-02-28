/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on March 19, 2021, 11:51 AM
 */

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
