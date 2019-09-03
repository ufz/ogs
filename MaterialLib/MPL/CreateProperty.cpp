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
    auto const property_type = config.getConfigParameter<std::string>("type");

    // If (and only if) the given property type is 'constant', a corresponding
    // value is needed.
    if (property_type == "Constant")
    {
        std::vector<double> const values =
            //! \ogs_file_param{properties__property__Constant__value}
            config.getConfigParameter<std::vector<double>>("value");

        switch (values.size())
        {
            case 1:
            {
                // scalar
                PropertyDataType property_value = values[0];
                return std::make_unique<Constant>(property_value);
            }
            case 2:
            {
                // Pair
                PropertyDataType property_value = Pair{values[0], values[1]};
                return std::make_unique<Constant>(property_value);
            }
            case 3:
            {
                // Vector
                PropertyDataType property_value =
                    Vector{values[0], values[1], values[2]};
                return std::make_unique<Constant>(property_value);
            }
            case 4:
            {
                // Tensor
                PropertyDataType property_value =
                    Tensor2d{values[0], values[1], values[2], values[3]};
                return std::make_unique<Constant>(property_value);
            }
            case 6:
            {
                // Symmetric Tensor - xx, yy, zz, xy, xz, yz
                PropertyDataType property_value =
                    SymmTensor{values[0], values[1], values[2],
                               values[3], values[4], values[5]};
                return std::make_unique<Constant>(property_value);
            }
            case 9:
            {
                // Tensor
                PropertyDataType property_value = Tensor{
                    values[0], values[1], values[2], values[3], values[4],
                    values[5], values[6], values[7], values[8]};
                return std::make_unique<Constant>(property_value);
            }

            default:
            {
                OGS_FATAL(
                    "Creation of a constant property with %i components is not "
                    "implemented.",
                    values.size());
            }
        }

        PropertyDataType property_value;
        return std::make_unique<Constant>(property_value);
    }
    // Properties can be medium, phase, or component properties.
    // Some of them require information about the respective material.
    // Thus, we pass a pointer to the material that requests the property.
    // In this method, this pointer is realized via typename MaterialType, which
    // replaces either Medium*, Phase*, or Component*.
    // Note that most property constructors (only those that request material
    // pointers) must be overloaded for any type of material.

    if (property_type == "Linear")
    {
        auto const reference_value =
            //! \ogs_file_param{properties__property__LinearProperty__reference_value}
            config.getConfigParameter<double>("reference_value");

        std::vector<MaterialPropertyLib::IndependentVariable> ivs;
        auto const& independent_variables_config =
            //! \ogs_file_param{properties__property__LinearProperty__independent_variable}
            config.getConfigSubtree("independent_variables");
        for (auto const& independent_variable_config :
             independent_variables_config.getConfigSubtreeList(
                 "independent_variable"))
        {
            auto const& variable_name =
                //! \ogs_file_param{properties__property__LinearProperty__independent_variable__variable_name}
                independent_variable_config.getConfigParameter<std::string>(
                    "variable_name");
            auto const reference_condition =
                //! \ogs_file_param{properties__property__LinearProperty__independent_variable__reference_condition}
                independent_variable_config.getConfigParameter<double>(
                    "reference_condition");
            auto const slope =
                //! \ogs_file_param{properties__property__LinearProperty__independent_variable__slope}
                independent_variable_config.getConfigParameter<double>("slope");

            MaterialPropertyLib::Variable ivt =
                MaterialPropertyLib::convertStringToVariable(variable_name);

            MaterialPropertyLib::IndependentVariable iv{
                ivt, reference_condition, slope};

            ivs.push_back(std::move(iv));
        }

        return std::make_unique<MaterialPropertyLib::LinearProperty>(
            reference_value, ivs);
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
            parameter_name, parameters, 1, nullptr);
        return std::make_unique<MaterialPropertyLib::ParameterProperty>(
            parameter);
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
