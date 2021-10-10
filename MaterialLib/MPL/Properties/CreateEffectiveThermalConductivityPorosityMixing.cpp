/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "BaseLib/ConfigTree.h"
#include "EffectiveThermalConductivityPorosityMixing.h"
#include "MaterialLib/MPL/Property.h"
#include "ParameterLib/CoordinateSystem.h"

namespace MaterialPropertyLib
{
std::unique_ptr<Property> createEffectiveThermalConductivityPorosityMixing(
    int const geometry_dimension,
    BaseLib::ConfigTree const& config,
    ParameterLib::CoordinateSystem const* const local_coordinate_system)
{
    //! \ogs_file_param{properties__property__type}
    config.checkConfigParameter("type",
                                "EffectiveThermalConductivityPorosityMixing");

    // Second access for storage.
    //! \ogs_file_param{properties__property__name}
    auto property_name = config.peekConfigParameter<std::string>("name");

    DBUG(
        "Create effective thermal_conductivity property from porosity mixing "
        "{:s}.",
        property_name);

    if (geometry_dimension == 1)
    {
        return std::make_unique<
            MaterialPropertyLib::EffectiveThermalConductivityPorosityMixing<1>>(
            std::move(property_name), local_coordinate_system);
    }

    if (geometry_dimension == 2)
    {
        return std::make_unique<
            MaterialPropertyLib::EffectiveThermalConductivityPorosityMixing<2>>(
            std::move(property_name), local_coordinate_system);
    }
    //! \ogs_file_param_special{properties__property__EffectiveThermalConductivityPorosityMixing}
    return std::make_unique<
        MaterialPropertyLib::EffectiveThermalConductivityPorosityMixing<3>>(
        std::move(property_name), local_coordinate_system);
}
}  // namespace MaterialPropertyLib
