// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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
