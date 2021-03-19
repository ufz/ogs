/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on March 4, 2021, 4:38 PM
 */

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

    // Here ogs_file_param is used just to create a documentation entry for this
    // property without any input parameter.
    //! \ogs_file_param{properties__property__WaterVapourDensity}
    auto property_name = config.peekConfigParameter<std::string>("name");

    return std::make_unique<WaterVapourDensity>(std::move(property_name));
}
}  // namespace MaterialPropertyLib
