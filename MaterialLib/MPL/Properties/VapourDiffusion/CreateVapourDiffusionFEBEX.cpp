/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on March 5, 2021, 4:51 PM
 */

#include "CreateVapourDiffusionFEBEX.h"

#include "BaseLib/ConfigTree.h"
#include "MaterialLib/MPL/Property.h"
#include "VapourDiffusionFEBEX.h"

namespace MaterialPropertyLib
{
std::unique_ptr<Property> createVapourDiffusionFEBEX(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{properties__property__type}
    config.checkConfigParameter("type", "VapourDiffusionFEBEX");
    DBUG("Create VapourDiffusionFEBEX medium property");

    // Second access for storage.
    //! \ogs_file_param{properties__property__name}
    auto property_name = config.peekConfigParameter<std::string>("name");

    auto const tortuosity =
        //! \ogs_file_param{properties__property__VapourDiffusionFEBEX__tortuosity}
        config.getConfigParameter<double>("tortuosity");

    return std::make_unique<VapourDiffusionFEBEX>(std::move(property_name),
                                                  tortuosity);
}
}  // namespace MaterialPropertyLib
