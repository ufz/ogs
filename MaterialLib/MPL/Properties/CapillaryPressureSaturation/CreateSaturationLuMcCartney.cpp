// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause
#include "BaseLib/ConfigTree.h"
#include "SaturationLuMcCartney.h"

namespace MaterialPropertyLib
{
std::unique_ptr<SaturationLuMcCartney> createSaturationLuMcCartney(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{properties__property__type}
    config.checkConfigParameter("type", "SaturationLuMcCartney");

    // Second access for storage.
    //! \ogs_file_param{properties__property__name}
    auto property_name = config.peekConfigParameter<std::string>("name");

    DBUG("Create SaturationLuMcCartney medium property {:s}.", property_name);

    auto const material =
        //! \ogs_file_param{properties__property__SaturationLuMcCartney__material}
        config.getConfigParameter<std::string>("material");

    return std::make_unique<SaturationLuMcCartney>(std::move(property_name),
                                                   material);
}
}  // namespace MaterialPropertyLib
