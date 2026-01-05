// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "CreateWaterVapourDensity.h"

#include "BaseLib/ConfigTree.h"
#include "MaterialLib/MPL/Property.h"
#include "WaterVapourDensity.h"

namespace MaterialPropertyLib
{
std::unique_ptr<Property> createWaterVapourDensity(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{properties__property__type}
    config.checkConfigParameter("type", "WaterVapourDensity");
    DBUG("Create WaterVapourDensity phase property");

    // Second access for storage.
    //! \ogs_file_param{properties__property__name}
    auto property_name = config.peekConfigParameter<std::string>("name");

    //! \ogs_file_param_special{properties__property__WaterVapourDensity}
    return std::make_unique<WaterVapourDensity>(std::move(property_name));
}
}  // namespace MaterialPropertyLib
