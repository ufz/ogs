/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on Feb 8, 2023, 4:38 PM
 */

#include "CreateWaterTemperatureIAPWSIF97Region1.h"

#include "BaseLib/ConfigTree.h"
#include "MaterialLib/MPL/Property.h"
#include "WaterTemperatureIAPWSIF97Region1.h"

namespace MaterialPropertyLib
{
std::unique_ptr<Property> createWaterTemperatureIAPWSIF97Region1(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{properties__property__type}
    config.checkConfigParameter("type", "WaterTemperatureIAPWSIF97Region1");
    DBUG("Create WaterTemperatureIAPWSIF97Region1 phase property");

    // Second access for storage.
    //! \ogs_file_param{properties__property__name}
    auto property_name = config.peekConfigParameter<std::string>("name");

    //! \ogs_file_param_special{properties__property__WaterTemperatureIAPWSIF97Region1}
    return std::make_unique<WaterTemperatureIAPWSIF97Region1>(
        std::move(property_name));
}
}  // namespace MaterialPropertyLib
