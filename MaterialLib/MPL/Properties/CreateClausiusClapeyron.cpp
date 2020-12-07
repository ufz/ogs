/**
 * \file
 * \author Norbert Grunwald
 * \date   Dec 07, 2020
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "BaseLib/ConfigTree.h"
#include "ClausiusClapeyron.h"

namespace MaterialPropertyLib
{
std::unique_ptr<ClausiusClapeyron> createClausiusClapeyron(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{properties__property__type}
    config.checkConfigParameter("type", "ClausiusClapeyron");

    //! \ogs_file_param{properties__property__name}
    auto property_name = config.peekConfigParameter<std::string>("name");

    DBUG("Create ClausiusClapeyron property {:s}.", property_name);

    auto const triple_temperature =
        //! \ogs_file_param{properties__property__ClausiusClapeyron__triple_temperature}
        config.getConfigParameter<double>("triple_temperature");

    auto const triple_pressure =
        //! \ogs_file_param{properties__property__ClausiusClapeyron__triple_pressure}
        config.getConfigParameter<double>("triple_pressure");

    auto const critical_temperature =
        //! \ogs_file_param{properties__property__ClausiusClapeyron__critical_temperature}
        config.getConfigParameter<double>("critical_temperature");

    auto const critical_pressure =
        //! \ogs_file_param{properties__property__ClausiusClapeyron__critical_pressure}
        config.getConfigParameter<double>("critical_pressure");

    auto const ref_temperature =
        //! \ogs_file_param{properties__property__ClausiusClapeyron__reference_temperature}
        config.getConfigParameter<double>("reference_temperature");

    auto const ref_pressure =
        //! \ogs_file_param{properties__property__ClausiusClapeyron__reference_pressure}
        config.getConfigParameter<double>("reference_pressure");

    //! \ogs_file_param_special{properties__property__ClausiusClapeyron}
    return std::make_unique<ClausiusClapeyron>(
        std::move(property_name), triple_temperature, triple_pressure,
        critical_temperature, critical_pressure, ref_temperature, ref_pressure);
}
}  // namespace MaterialPropertyLib
