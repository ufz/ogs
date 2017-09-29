/**
 * \author Norbert Grunwald
 * \date   Sep 7, 2017
 * \brief
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "mpComponent.h"
#include "Components/components.h"
#include "Properties/properties.h"

namespace MaterialPropertyLib
{
/// Default constructor of Component. This constructor is used
/// when the component is not specified via the 'name'-tag.
Component::Component()
{
    // DBUG ("      Unspecified component constructed.");
    // Default properties are set to zero.
    createDefaultProperties();

    // Some properties can be initialized by other default
    // properties:
    _properties[name] = std::make_unique<Constant>("no_name");
}
/**
 * \brief The method reads the 'properties' tag in the prj-file and creates
 * component properties accordingly.
 * \details First, a new property iy created based on the specified property
 * type. Then, the property name is evaluated and the property is copied
 * into the properties array.
 */
void Component::createProperties(BaseLib::ConfigTree const& config)
{
    for (auto property_config : config.getConfigSubtreeList("property"))
    {
        // Parsing the property name:
        auto const property_name =
            property_config.getConfigParameter<std::string>("name");
        // Create a new property based on the configuration subtree:
        auto property = newProperty(property_config, this);

        // Insert the new property at the right position into the components
        // private PropertyArray:
        _properties[convertStringToProperty(property_name)] =
            std::move(property);
    }
}

/**
 * This method initializes the component property array with
 * constant functions of value zero.
 */
void Component::createDefaultProperties(void)
{
    for (std::size_t i = 0; i < number_of_property_enums; ++i)
        _properties[i] = std::make_unique<Constant>(0);
}

Property* Component::property(PropertyEnum const& p)
{
    return _properties[p].get();
}

void Component::resetPropertyUpdateStatus(void)
{
    for (std::size_t p = 0; p < number_of_property_enums; ++p)
        _properties[p]->isUpdated(false);
}
/**
 * \brief This function creates a new component based on the (optional)
 * component name that is given in the prj-file.
 * \details The method evaluates the string in the 'name'-object and
 * calls the constructors of the derived component classes (if found)
 * or that of the base class (if no name is specified).
 */
std::unique_ptr<Component> newComponent(
    boost::optional<std::string> const& name)
{
    // Check whether a name is given or not
    if (!name)
        // If there is no component name given (which is common), the
        // method creates a custom component, where all properties have
        // to be specified manually.
        return std::make_unique<Component>();

    std::string component_name = static_cast<std::string>(name.get());
    // If a name is given, it must conform with one of the
    // derived component names in the following list:
    if (boost::iequals(component_name, "water"))  // "string" == "StRiNg"
    {
        return std::make_unique<Water>();
    }
    if (boost::iequals(component_name, "sodiumchloride"))
    {
        return std::make_unique<SodiumChloride>();
    }
    if (boost::iequals(component_name, "carbondioxide"))
    {
        return std::make_unique<CarbonDioxide>();
    }
    // If the given name does not appear in the list, this
    // throws a fatal error.
    OGS_FATAL("The specified component name \"%s\" was not recognized",
              component_name.c_str());
    return nullptr;  // to avoid the 'no return' warning.
}
}  // MaterialPropertyLib
