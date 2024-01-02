/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

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

    return std::make_unique<MaterialPropertyLib::TemperatureDependentFraction>(
        std::move(property_name), k, T_c);
}
}  // namespace MaterialPropertyLib
