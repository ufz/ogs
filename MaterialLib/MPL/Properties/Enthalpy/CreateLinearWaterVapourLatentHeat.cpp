/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on March 16, 2021, 10:13 AM
 */

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
