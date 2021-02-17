/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on February 17, 2021, 3:47 PM
 */

#include "CreateSoilThermalConductivitySomerton.h"

#include "BaseLib/ConfigTree.h"
#include "MaterialLib/MPL/Property.h"
#include "SoilThermalConductivitySomerton.h"

namespace MaterialPropertyLib
{
std::unique_ptr<Property> createSoilThermalConductivitySomerton(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{properties__property__type}
    config.checkConfigParameter("type", "SoilThermalConductivitySomerton");

    //! \ogs_file_param{properties__property__name}
    auto property_name = config.peekConfigParameter<std::string>("name");

    DBUG("Create SoilThermalConductivitySomerton medium property");

    auto const dry_thermal_conductivity =
        //! \ogs_file_param{properties__property__SoilThermalConductivitySomerton__dry_thermal_conductivity}
        config.getConfigParameter<double>("dry_thermal_conductivity");

    auto const wet_thermal_conductivity =
        //! \ogs_file_param{properties__property__SoilThermalConductivitySomerton__wet_thermal_conductivity}
        config.getConfigParameter<double>("wet_thermal_conductivity");

    return std::make_unique<SoilThermalConductivitySomerton>(
        property_name, dry_thermal_conductivity, wet_thermal_conductivity);
}
}  // namespace MaterialPropertyLib
