/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on Feb 8, 2023, 12:31 PM
 */

#include "CreateWaterVapourEnthalpyIAPWSIF97Region4.h"

#include "BaseLib/ConfigTree.h"
#include "MaterialLib/MPL/Property.h"
#include "WaterVapourEnthalpyIAPWSIF97Region4.h"

namespace MaterialPropertyLib
{
std::unique_ptr<Property> createWaterVapourEnthalpyIAPWSIF97Region4(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{properties__property__type}
    config.checkConfigParameter("type", "WaterVapourEnthalpyIAPWSIF97Region4");
    DBUG("Create WaterVapourEnthalpyIAPWSIF97Region4 phase property");

    // Second access for storage.
    //! \ogs_file_param{properties__property__name}
    auto property_name = config.peekConfigParameter<std::string>("name");

    //! \ogs_file_param_special{properties__property__WaterVapourEnthalpyIAPWSIF97Region4}
    return std::make_unique<WaterVapourEnthalpyIAPWSIF97Region4>(
        std::move(property_name));
}
}  // namespace MaterialPropertyLib
