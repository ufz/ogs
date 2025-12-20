// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "CreateLinearWaterVapourLatentHeat.h"

#include "BaseLib/ConfigTree.h"
#include "LinearWaterVapourLatentHeat.h"
#include "MaterialLib/MPL/Property.h"

namespace MaterialPropertyLib
{
std::unique_ptr<Property> createLinearWaterVapourLatentHeat(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{properties__property__type}
    config.checkConfigParameter("type", "LinearWaterVapourLatentHeat");
    DBUG("Create LinearWaterVapourLatentHeat phase property");

    // Second access for storage.
    //! \ogs_file_param{properties__property__name}
    auto property_name = config.peekConfigParameter<std::string>("name");

    //! \ogs_file_param_special{properties__property__LinearWaterVapourLatentHeat}
    return std::make_unique<LinearWaterVapourLatentHeat>(
        std::move(property_name));
}
}  // namespace MaterialPropertyLib
