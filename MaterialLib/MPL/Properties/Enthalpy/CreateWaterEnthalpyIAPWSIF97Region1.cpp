// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "CreateWaterEnthalpyIAPWSIF97Region1.h"

#include "BaseLib/ConfigTree.h"
#include "MaterialLib/MPL/Property.h"
#include "WaterEnthalpyIAPWSIF97Region1.h"

namespace MaterialPropertyLib
{
std::unique_ptr<Property> createWaterEnthalpyIAPWSIF97Region1(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{properties__property__type}
    config.checkConfigParameter("type", "WaterEnthalpyIAPWSIF97Region1");
    DBUG("Create WaterEnthalpyIAPWSIF97Region1 phase property");

    // Second access for storage.
    //! \ogs_file_param{properties__property__name}
    auto property_name = config.peekConfigParameter<std::string>("name");

    //! \ogs_file_param_special{properties__property__WaterEnthalpyIAPWSIF97Region1}
    return std::make_unique<WaterEnthalpyIAPWSIF97Region1>(
        std::move(property_name));
}
}  // namespace MaterialPropertyLib
