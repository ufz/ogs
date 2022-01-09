/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "BaseLib/ConfigTree.h"
#include "ParameterLib/CoordinateSystem.h"
#include "SaturationDependentSwelling.h"

namespace MaterialPropertyLib
{
std::unique_ptr<SaturationDependentSwelling> createSaturationDependentSwelling(
    BaseLib::ConfigTree const& config,
    ParameterLib::CoordinateSystem const* const local_coordinate_system)
{
    //! \ogs_file_param{properties__property__type}
    config.checkConfigParameter("type", "SaturationDependentSwelling");

    // Second access for storage.
    //! \ogs_file_param{properties__property__name}
    auto property_name = config.peekConfigParameter<std::string>("name");

    DBUG("Create SaturationDependentSwelling solid phase property {:s}.",
         property_name);

    auto const swelling_pressures =
        //! \ogs_file_param{properties__property__SaturationDependentSwelling__swelling_pressures}
        config.getConfigParameter<std::vector<double>>("swelling_pressures");

    if (swelling_pressures.size() != 3)
    {
        OGS_FATAL(
            "The number of swelling pressures must be three, but {:d} were "
            "given.",
            swelling_pressures.size());
    }

    auto const exponents =
        //! \ogs_file_param{properties__property__SaturationDependentSwelling__exponents}
        config.getConfigParameter<std::vector<double>>("exponents");

    if (exponents.size() != 3)
    {
        OGS_FATAL("The number of exponents must be three, but {:d} were given.",
                  exponents.size());
    }

    if (swelling_pressures.size() != exponents.size())
    {
        OGS_FATAL(
            "The number of swelling pressures and exponents must be equal, but "
            "they are {:d} and {:d}, respectively.",
            swelling_pressures.size(), exponents.size());
    }

    auto const lower_saturation_limit =
        //! \ogs_file_param{properties__property__SaturationDependentSwelling__lower_saturation_limit}
        config.getConfigParameter<double>("lower_saturation_limit");

    auto const upper_saturation_limit =
        //! \ogs_file_param{properties__property__SaturationDependentSwelling__upper_saturation_limit}
        config.getConfigParameter<double>("upper_saturation_limit");

    return std::make_unique<SaturationDependentSwelling>(
        std::move(property_name),
        std::array<double, 3>{swelling_pressures[0], swelling_pressures[1],
                              swelling_pressures[2]},
        std::array<double, 3>{exponents[0], exponents[1], exponents[2]},
        lower_saturation_limit, upper_saturation_limit,
        local_coordinate_system);
}
}  // namespace MaterialPropertyLib
