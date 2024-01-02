/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on February 17, 2021, 3:47 PM
 */

#include "CreateSaturationWeightedThermalConductivity.h"

#include <boost/mp11.hpp>
#include <string>

#include "BaseLib/Algorithm.h"
#include "BaseLib/ConfigTree.h"
#include "MaterialLib/MPL/Property.h"
#include "ParameterLib/Parameter.h"
#include "ParameterLib/Utils.h"
#include "SaturationWeightedThermalConductivity.h"

namespace
{
template <MaterialPropertyLib::MeanType MeanType, int Dim>
std::unique_ptr<MaterialPropertyLib::Property>
createSaturationWeightedThermalConductivity(
    std::string name,
    ParameterLib::Parameter<double> const& dry_thermal_conductivity,
    ParameterLib::Parameter<double> const& wet_thermal_conductivity)
{
    return std::make_unique<
        MaterialPropertyLib::SaturationWeightedThermalConductivity<MeanType,
                                                                   Dim>>(
        std::move(name), dry_thermal_conductivity, wet_thermal_conductivity);
}
}  // namespace

namespace MaterialPropertyLib
{
std::unique_ptr<Property> createSaturationWeightedThermalConductivity(
    int const geometry_dimension,
    BaseLib::ConfigTree const& config,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters)
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

    std::map<
        std::pair<MeanType, int>,
        std::unique_ptr<Property> (*)(
            std::string /*name*/,
            ParameterLib::Parameter<double> const& /*dry_thermal_conductivity*/,
            ParameterLib::Parameter<
                double> const& /*wet_thermal_conductivity*/)>
        map_dim_and_mean_to_creator;

    // initialize the map
    {
        using namespace boost::mp11;
        using Dims = mp_list<mp_int<1>, mp_int<2>, mp_int<3>>;
        using Means = mp_list<
            std::integral_constant<MeanType, MeanType::ARITHMETIC_LINEAR>,
            std::integral_constant<MeanType, MeanType::ARITHMETIC_SQUAREROOT>,
            std::integral_constant<MeanType, MeanType::GEOMETRIC>>;
        using DimsAndMeanTypes =
            mp_product<mp_list, Dims,
                       Means>;  // Cartesian product of Dims and Means.

        mp_for_each<DimsAndMeanTypes>(
            [&map_dim_and_mean_to_creator]<typename Dim, typename Mean>(
                mp_list<Dim, Mean>)
            {
                map_dim_and_mean_to_creator.emplace(
                    std::pair{Mean::value, Dim::value},
                    &::createSaturationWeightedThermalConductivity<Mean::value,
                                                                   Dim::value>);
            });
    }

    auto const it = map_dim_and_mean_to_creator.find(
        std::pair{mean_type, geometry_dimension});

    if (it == map_dim_and_mean_to_creator.end())
    {
        OGS_FATAL(
            "Cannot create a SaturationWeightedThermalConductivity model for "
            "dimension {} and mean type {}",
            geometry_dimension, mean_type_str);
    }

    auto* creator = it->second;
    return creator(std::move(property_name),
                   dry_thermal_conductivity,
                   wet_thermal_conductivity);
}
}  // namespace MaterialPropertyLib
