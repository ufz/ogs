/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on March 16, 2021, 10:13 AM
 */

#include "CreateWaterVapourLatentHeat.h"

#include "BaseLib/ConfigTree.h"
#include "MaterialLib/MPL/Property.h"
#include "WaterVapourLatentHeat.h"

namespace MaterialPropertyLib
{
std::unique_ptr<Property> createWaterVapourLatentHeat(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{properties__property__type}
    config.checkConfigParameter("type", "WaterVapourLatentHeat");
    DBUG("Create WaterVapourLatentHeat medium property");

    // Second access for storage.
    //! \ogs_file_param{properties__property__name}
    auto property_name = config.peekConfigParameter<std::string>("name");

    //! \ogs_file_param_special{properties__property__WaterVapourLatentHeat}
    return std::make_unique<WaterVapourLatentHeat>(std::move(property_name));
}
}  // namespace MaterialPropertyLib
