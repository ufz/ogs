// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "BaseLib/ConfigTree.h"
#include "TemperatureDependentFraction.h"

namespace MaterialPropertyLib
{
std::unique_ptr<TemperatureDependentFraction>
createTemperatureDependentFraction(BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{properties__property__type}
    config.checkConfigParameter("type", "TemperatureDependentFraction");

    // Second access for storage.
    //! \ogs_file_param{properties__property__name}
    auto property_name = config.peekConfigParameter<std::string>("name");

    DBUG("Create temperature dependent fraction property {:s}.", property_name);

    auto const k =
        //! \ogs_file_param{properties__property__TemperatureDependentFraction__steepness}
        config.getConfigParameter<double>("steepness");

    auto const T_c =
        //! \ogs_file_param{properties__property__TemperatureDependentFraction__characteristic_temperature}
        config.getConfigParameter<double>("characteristic_temperature");

    auto const S_r =
        //! \ogs_file_param{properties__property__TemperatureDependentFraction__residual_saturation}
        config.getConfigParameter<double>("residual_saturation", 0.0);

    if (S_r < 0.0 || S_r > 1.0)
    {
        OGS_FATAL(
            "In the temperature dependent fraction property setting, "
            "the residual_saturation value {} is out of [0, 1]-range.",
            S_r);
    }

    return std::make_unique<MaterialPropertyLib::TemperatureDependentFraction>(
        std::move(property_name), k, T_c, S_r);
}
}  // namespace MaterialPropertyLib
