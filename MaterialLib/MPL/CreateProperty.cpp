/**
 * \file
 * \author Norbert Grunwald
 * \date   Sep 7, 2017
 *
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateProperty.h"

#include <string>
#include <vector>
#include "BaseLib/ConfigTree.h"
#include "ParameterLib/Parameter.h"
#include "ParameterLib/Utils.h"

#include "Properties/CreateProperties.h"
#include "Properties/Properties.h"

#include "Component.h"
#include "Medium.h"
#include "Phase.h"

namespace
{
std::unique_ptr<MaterialPropertyLib::Property> createProperty(
    BaseLib::ConfigTree const& config,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters)
{
    using namespace MaterialPropertyLib;
    // Parsing the property type:
    //! \ogs_file_param{properties__property__type}
    auto const property_type = config.peekConfigParameter<std::string>("type");

    if (property_type == "Constant")
    {
        return createConstant(config);
    }
    if (property_type == "Linear")
    {
        return createLinearProperty(config);
    }

    if (property_type == "Exponential")
    {
        auto const reference_value =
            //! \ogs_file_param{properties__property__ExponentialProperty__reference_value}
            config.getConfigParameter<double>("reference_value");

        auto const& exponent_data_config =
            //! \ogs_file_param{properties__property__ExponentialProperty__exponent}
            config.getConfigSubtree("exponent");

        auto const& variable_name =
            //! \ogs_file_param{properties__property__ExponentialProperty__exponent__variable_name}
            exponent_data_config.getConfigParameter<std::string>(
                "variable_name");
        auto const reference_condition =
            //! \ogs_file_param{properties__property__ExponentialProperty__exponent__reference_condition}
            exponent_data_config.getConfigParameter<double>(
                "reference_condition");
        auto const factor =
            //! \ogs_file_param{properties__property__ExponentialProperty__exponent__factor}
            exponent_data_config.getConfigParameter<double>("factor");

        MaterialPropertyLib::Variable exp_data_type =
            MaterialPropertyLib::convertStringToVariable(variable_name);

        MaterialPropertyLib::ExponentData const exp_data{
            exp_data_type, reference_condition, factor};

        return std::make_unique<MaterialPropertyLib::ExponentialProperty>(
            reference_value, exp_data);
    }

    if (property_type == "Parameter")
    {
        std::string const& parameter_name =
            //! \ogs_file_param{properties__property__Parameter__parameter_name}
            config.getConfigParameter<std::string>("parameter_name");
        auto const& parameter = ParameterLib::findParameter<double>(
            parameter_name, parameters, 0, nullptr);
        return std::make_unique<MaterialPropertyLib::ParameterProperty>(
            parameter);
    }

    if (boost::iequals(property_type, "IdealGasLaw"))
    {
        return createIdealGasLaw(config);
    }
    /* TODO Additional properties go here, for example:
    if (boost::iequals(property_type, "BilinearTemperaturePressure"))
    {
        return createBilinearTemperaturePressure(config, material_type);
    }
    */

    // If none of the above property types are found, OGS throws an error.
    OGS_FATAL("The specified component property type '%s' was not recognized",
              property_type.c_str());
}
}  // namespace

namespace MaterialPropertyLib
{
std::unique_ptr<PropertyArray> createProperties(
    boost::optional<BaseLib::ConfigTree> const& config,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters)
{
    if (!config)
    {
        return nullptr;
    }

    //! \ogs_file_param{properties__property}
    auto const& property_configs = config->getConfigSubtreeList("property");
    if (property_configs.empty())
    {
        return nullptr;
    }

    auto properties = std::make_unique<PropertyArray>();

    for (auto property_config : property_configs)
    {
        // Parsing the property name:
        auto const property_name =
            //! \ogs_file_param{properties__property__name}
            property_config.getConfigParameter<std::string>("name");
        // Create a new property based on the configuration subtree:
        auto property = createProperty(property_config, parameters);

        // Insert the new property at the right position into the components
        // private PropertyArray:
        (*properties)[convertStringToProperty(property_name)] =
            std::move(property);
    }
    return properties;
}

}  // namespace MaterialPropertyLib
