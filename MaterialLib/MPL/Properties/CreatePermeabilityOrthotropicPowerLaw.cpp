/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "BaseLib/ConfigTree.h"
#include "ParameterLib/CoordinateSystem.h"
#include "PermeabilityOrthotropicPowerLaw.h"

namespace MaterialPropertyLib
{
std::unique_ptr<Property> createPermeabilityOrthotropicPowerLaw(
    BaseLib::ConfigTree const& config,
    ParameterLib::CoordinateSystem const* const local_coordinate_system)
{
    //! \ogs_file_param{properties__property__type}
    config.checkConfigParameter("type", "PermeabilityOrthotropicPowerLaw");

    // Second access for storage.
    //! \ogs_file_param{properties__property__name}
    auto property_name = config.peekConfigParameter<std::string>("name");

    DBUG("Create PermeabilityOrthotropicPowerLaw solid phase property {:s}.",
         property_name);

    auto const intrinsic_permeabilities =
        //! \ogs_file_param{properties__property__PermeabilityOrthotropicPowerLaw__intrinsic_permeabilities}
        config.getConfigParameter<std::vector<double>>(
            "intrinsic_permeabilities");

    if (!((intrinsic_permeabilities.size() == 3) ||
          (intrinsic_permeabilities.size() == 2)))
    {
        OGS_FATAL(
            "The number of intrinsic permeabilities must be two or three, but "
            "{:d} were given.",
            intrinsic_permeabilities.size());
    }

    auto const exponents =
        //! \ogs_file_param{properties__property__PermeabilityOrthotropicPowerLaw__exponents}
        config.getConfigParameter<std::vector<double>>("exponents");

    if (exponents.size() != 3 && exponents.size() != 2)
    {
        OGS_FATAL(
            "The number of exponents must be two or three, but {:d} were "
            "given.",
            exponents.size());
    }

    if (intrinsic_permeabilities.size() != exponents.size())
    {
        OGS_FATAL(
            "The number of intrinsic permeabilities and exponents must be "
            "equal, but they are {:d} and {:d}, respectively.",
            intrinsic_permeabilities.size(), exponents.size());
    }

    if (exponents.size() == 2)
    {
        return std::make_unique<PermeabilityOrthotropicPowerLaw<2>>(
            std::move(property_name),
            std::array<double, 2>{intrinsic_permeabilities[0],
                                  intrinsic_permeabilities[1]},
            std::array<double, 2>{exponents[0], exponents[1]},
            local_coordinate_system);
    }
    if (exponents.size() == 3)
    {
        return std::make_unique<PermeabilityOrthotropicPowerLaw<3>>(
            std::move(property_name),
            std::array<double, 3>{intrinsic_permeabilities[0],
                                  intrinsic_permeabilities[1],
                                  intrinsic_permeabilities[2]},
            std::array<double, 3>{exponents[0], exponents[1], exponents[2]},
            local_coordinate_system);
    }
    OGS_FATAL(
        "Could not create PermeabilityOrthotropicPowerLaw material model.");
}
}  // namespace MaterialPropertyLib
