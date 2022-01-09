/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on February 17, 2021, 3:47 PM
 */

#include "CreateSoilThermalConductivitySomerton.h"

#include <string>

#include "BaseLib/ConfigTree.h"
#include "MaterialLib/MPL/Property.h"
#include "ParameterLib/CoordinateSystem.h"
#include "ParameterLib/Parameter.h"
#include "ParameterLib/Utils.h"
#include "SoilThermalConductivitySomerton.h"

namespace MaterialPropertyLib
{
std::unique_ptr<Property> createSoilThermalConductivitySomerton(
    int const geometry_dimension,
    BaseLib::ConfigTree const& config,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    ParameterLib::CoordinateSystem const* const local_coordinate_system)
{
    //! \ogs_file_param{properties__property__type}
    config.checkConfigParameter("type", "SoilThermalConductivitySomerton");

    //! \ogs_file_param{properties__property__name}
    auto property_name = config.peekConfigParameter<std::string>("name");

    DBUG("Create SoilThermalConductivitySomerton medium property");

    std::string const& dry_thermal_conductivity_parameter_name =
        //! \ogs_file_param{properties__property__SoilThermalConductivitySomerton__dry_thermal_conductivity}
        config.getConfigParameter<std::string>("dry_thermal_conductivity");
    auto const& dry_thermal_conductivity = ParameterLib::findParameter<double>(
        dry_thermal_conductivity_parameter_name, parameters, 0, nullptr);

    std::string const& wet_thermal_conductivity_parameter_name =
        //! \ogs_file_param{properties__property__SoilThermalConductivitySomerton__wet_thermal_conductivity}
        config.getConfigParameter<std::string>("wet_thermal_conductivity");
    auto const& wet_thermal_conductivity = ParameterLib::findParameter<double>(
        wet_thermal_conductivity_parameter_name, parameters, 0, nullptr);

    if (geometry_dimension == 1)
    {
        return std::make_unique<SoilThermalConductivitySomerton<1>>(
            std::move(property_name),
            dry_thermal_conductivity,
            wet_thermal_conductivity,
            local_coordinate_system);
    }

    if (geometry_dimension == 2)
    {
        return std::make_unique<SoilThermalConductivitySomerton<2>>(
            std::move(property_name),
            dry_thermal_conductivity,
            wet_thermal_conductivity,
            local_coordinate_system);
    }

    return std::make_unique<SoilThermalConductivitySomerton<3>>(
        std::move(property_name),
        dry_thermal_conductivity,
        wet_thermal_conductivity,
        local_coordinate_system);
}
}  // namespace MaterialPropertyLib
