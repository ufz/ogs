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

#include "CreateSaturationWeightedThermalConductivity.h"

#include <string>

#include "BaseLib/Algorithm.h"
#include "BaseLib/ConfigTree.h"
#include "MaterialLib/MPL/Property.h"
#include "ParameterLib/CoordinateSystem.h"
#include "ParameterLib/Parameter.h"
#include "ParameterLib/Utils.h"
#include "SaturationWeightedThermalConductivity.h"

namespace MaterialPropertyLib
{
std::unique_ptr<Property> createSaturationWeightedThermalConductivity(
    int const geometry_dimension,
    BaseLib::ConfigTree const& config,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    ParameterLib::CoordinateSystem const* const local_coordinate_system)
{
    //! \ogs_file_param{properties__property__type}
    config.checkConfigParameter("type",
                                "SaturationWeightedThermalConductivity");

    //! \ogs_file_param{properties__property__name}
    auto property_name = config.peekConfigParameter<std::string>("name");

    DBUG("Create SaturationWeightedThermalConductivity medium property");

    std::string const& dry_thermal_conductivity_parameter_name =
        //! \ogs_file_param{properties__property__SaturationWeightedThermalConductivity__dry_thermal_conductivity}
        config.getConfigParameter<std::string>("dry_thermal_conductivity");
    auto const& dry_thermal_conductivity = ParameterLib::findParameter<double>(
        dry_thermal_conductivity_parameter_name, parameters, 0, nullptr);

    std::string const& wet_thermal_conductivity_parameter_name =
        //! \ogs_file_param{properties__property__SaturationWeightedThermalConductivity__wet_thermal_conductivity}
        config.getConfigParameter<std::string>("wet_thermal_conductivity");
    auto const& wet_thermal_conductivity = ParameterLib::findParameter<double>(
        wet_thermal_conductivity_parameter_name, parameters, 0, nullptr);

    std::string const& mean_type_str =
        //! \ogs_file_param{properties__property__SaturationWeightedThermalConductivity__mean_type}
        config.getConfigParameter<std::string>("mean_type");

    const std::map<std::string, MeanType> mean_type_map{
        {"arithmetic_linear", MeanType::ARITHMETIC_LINEAR},
        {"arithmetic_squareroot", MeanType::ARITHMETIC_SQUAREROOT},
        {"geometric", MeanType::GEOMETRIC}};
    MeanType const& mean_type = BaseLib::getOrError(
        mean_type_map, mean_type_str,
        "Specified mean type for the thermal conductivity could not be found.");

    switch (geometry_dimension)
    {
        case 1:
        {
            switch (mean_type)
            {
                case MeanType::ARITHMETIC_LINEAR:
                {
                    return std::make_unique<SaturationWeightedThermalConductivity<
                        MeanType::ARITHMETIC_LINEAR, 1>>(
                            std::move(property_name),
                            dry_thermal_conductivity,
                            wet_thermal_conductivity,
                            local_coordinate_system);
                }
                case MeanType::ARITHMETIC_SQUAREROOT:
                {
                    return std::make_unique<SaturationWeightedThermalConductivity<
                        MeanType::ARITHMETIC_SQUAREROOT, 1>>(
                            std::move(property_name),
                            dry_thermal_conductivity,
                            wet_thermal_conductivity,
                            local_coordinate_system);
                }
                case MeanType::GEOMETRIC:
                {
                    return std::make_unique<SaturationWeightedThermalConductivity<
                        MeanType::GEOMETRIC, 1>>(
                            std::move(property_name),
                            dry_thermal_conductivity,
                            wet_thermal_conductivity,
                            local_coordinate_system);
                }
                default:
                {
                    OGS_FATAL("The requested MeanType is not implemented.");
                }
            }
        }
        case 2:
        {
            switch (mean_type)
            {
                case MeanType::ARITHMETIC_LINEAR:
                {
                    return std::make_unique<SaturationWeightedThermalConductivity<
                        MeanType::ARITHMETIC_LINEAR, 2>>(
                            std::move(property_name),
                            dry_thermal_conductivity,
                            wet_thermal_conductivity,
                            local_coordinate_system);
                }
                case MeanType::ARITHMETIC_SQUAREROOT:
                {
                    return std::make_unique<SaturationWeightedThermalConductivity<
                        MeanType::ARITHMETIC_SQUAREROOT, 2>>(
                            std::move(property_name),
                            dry_thermal_conductivity,
                            wet_thermal_conductivity,
                            local_coordinate_system);
                }
                case MeanType::GEOMETRIC:
                {
                    return std::make_unique<SaturationWeightedThermalConductivity<
                        MeanType::GEOMETRIC, 2>>(
                            std::move(property_name),
                            dry_thermal_conductivity,
                            wet_thermal_conductivity,
                            local_coordinate_system);
                }
                default:
                {
                    OGS_FATAL("The requested MeanType is not implemented.");
                }
            }

        }
        case 3:
        {
            switch (mean_type)
            {
                case MeanType::ARITHMETIC_LINEAR:
                {
                    return std::make_unique<SaturationWeightedThermalConductivity<
                        MeanType::ARITHMETIC_LINEAR, 3>>(
                            std::move(property_name),
                            dry_thermal_conductivity,
                            wet_thermal_conductivity,
                            local_coordinate_system);
                }
                case MeanType::ARITHMETIC_SQUAREROOT:
                {
                    return std::make_unique<SaturationWeightedThermalConductivity<
                        MeanType::ARITHMETIC_SQUAREROOT, 3>>(
                            std::move(property_name),
                            dry_thermal_conductivity,
                            wet_thermal_conductivity,
                            local_coordinate_system);
                }
                case MeanType::GEOMETRIC:
                {
                    return std::make_unique<SaturationWeightedThermalConductivity<
                        MeanType::GEOMETRIC, 3>>(
                            std::move(property_name),
                            dry_thermal_conductivity,
                            wet_thermal_conductivity,
                            local_coordinate_system);
                }
                default:
                {
                    OGS_FATAL("The requested MeanType is not implemented.");
                }
            }
        }
        default:
        {
            OGS_FATAL("Dimension doesn't exist.");
        }
    }

}
}  // namespace MaterialPropertyLib
