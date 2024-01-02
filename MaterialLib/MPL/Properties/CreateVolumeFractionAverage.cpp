/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "BaseLib/ConfigTree.h"
#include "VolumeFractionAverage.h"

namespace MaterialPropertyLib
{
std::unique_ptr<VolumeFractionAverage> createVolumeFractionAverage(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{properties__property__type}
    config.checkConfigParameter("type", "VolumeFractionAverage");

    // Second access for storage.
    //! \ogs_file_param{properties__property__name}
    auto property_name = config.peekConfigParameter<std::string>("name");

    //! \ogs_file_param_special{properties__property__VolumeFractionAverage}
    DBUG("Create volume fraction average {:s}.", property_name);

    // no input parameters required here (taken from phase properties)

    return std::make_unique<MaterialPropertyLib::VolumeFractionAverage>(
        std::move(property_name));
}
}  // namespace MaterialPropertyLib
