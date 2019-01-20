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

#include "Properties/Properties.h"

#include "Component.h"
#include "Medium.h"
#include "Phase.h"

namespace
{
std::unique_ptr<MaterialPropertyLib::Property> createProperty(
    BaseLib::ConfigTree const& config)
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
    boost::optional<BaseLib::ConfigTree> const& config)
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
        auto property = createProperty(property_config);

        // Insert the new property at the right position into the components
        // private PropertyArray:
        (*properties)[convertStringToProperty(property_name)] =
            std::move(property);
    }
    return properties;
}

}  // namespace MaterialPropertyLib
