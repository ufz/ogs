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

#include "mpProperty.h"
#include "Properties/properties.h"
#include "mpComponent.h"
#include "mpMedium.h"
#include "mpPhase.h"

namespace MaterialPropertyLib
{
/// The base class constructor may remain empty since the base class is
/// almost abstract.
Property::Property() : _isUpdated(0){};

PropertyDataType Property::value() const
{
    return _value;
}
/// The default implementation of this method only returns the
/// property value without altering it.
PropertyDataType Property::value(VariableArray const&)
{
    return _value;
}

void Property::isUpdated(bool status)
{
    _isUpdated = status;
}

bool Property::isUpdated()
{
    return _isUpdated;
}

void Property::notImplemented(std::string property, std::string material)
{
    OGS_FATAL(
        "The property \'%s\' is not "
        "available on the \'%s\' scale",
        property, material);
}

std::unique_ptr<Property> newProperty(BaseLib::ConfigTree const& config,
                                      Medium* m)
{
    return selectProperty(config, m);
}

std::unique_ptr<Property> newProperty(BaseLib::ConfigTree const& config,
                                      Phase* p)
{
    return selectProperty(config, p);
}

std::unique_ptr<Property> newProperty(BaseLib::ConfigTree const& config,
                                      Component* c)
{
    return selectProperty(config, c);
}

/**
 * This method creates a new Property. It uses the information stored
 * inside the configuration tree to select the property that is specified
 * by the user. It further passes a pointer to the material that requests
 * the property.
*/
template <typename MaterialType>
std::unique_ptr<Property> selectProperty(BaseLib::ConfigTree const& config,
                                         MaterialType M)
{
    // Parsing the property type:
    auto const property_type = config.getConfigParameter<std::string>("type");

    // If (and only if) the given property type is 'constant', a
    // corresponding value is needed.
    if (boost::iequals(property_type, "constant"))
    {
        /// \todo: Creating constant properties from prj-file is
        /// currently restricted to scalar values. The Constant
        /// constructor, however, can handle any datatype defined
        /// by PropertyDataType. This could be enhanced in order
        /// to define vectors or even tensors as constant properties.
        auto const property_value = config.getConfigParameter<double>("value");
        return std::make_unique<Constant>(Constant(property_value));
    }
    /// Properties can be medium, phase, or component properties.
    /// Some of them require information about the respective material.
    /// Thus, we pass a pointer to the material that requests the property.
    /// In this method, this pointer is realized via typename MaterialType,
    /// which replaces either Medium*, Phase*, or Component*.
    /// Note that most property constructors (only those that request material
    /// pointers) must be overloaded for any type of material.
    if (boost::iequals(property_type, "LinearTemperature"))
    {
        return std::make_unique<LinearTemperature>(M);
    }
    if (boost::iequals(property_type, "LinearEpsilon"))
    {
        DBUG("TODO: Implementation of %s property!!", property_type.c_str());
        return std::make_unique<Constant>(0);
    }
    if (boost::iequals(property_type, "AverageMolarMass"))
    {
        DBUG("TODO: Implementation of %s property!!", property_type.c_str());
        return std::make_unique<AverageMoleFraction>(M);
    }
    if (boost::iequals(property_type, "AverageVolumeFraction"))
    {
        DBUG("TODO: Implementation of %s property!!", property_type.c_str());
        return std::make_unique<Constant>(0);
    }
    if (boost::iequals(property_type, "Duan_2012"))
    {
        DBUG("TODO: Implementation of %s property!!", property_type.c_str());
        return std::make_unique<Constant>(0);
    }
    if (boost::iequals(property_type, "IAPWS_2008"))
    {
        DBUG("TODO: Implementation of %s property!!", property_type.c_str());
        return std::make_unique<Constant>(0);
    }
    if (boost::iequals(property_type, "Islam_Carlson_2012"))
    {
        DBUG("TODO: Implementation of %s property!!", property_type.c_str());
        return std::make_unique<Constant>(0);
    }
    if (boost::iequals(property_type, "Fenghour_1998"))
    {
        DBUG("TODO: Implementation of %s property!!", property_type.c_str());
        return std::make_unique<Constant>(0);
    }
    if (boost::iequals(property_type, "Mualem_1978"))
    {
        DBUG("TODO: Implementation of %s property!!", property_type.c_str());
        return std::make_unique<Constant>(0);
    }
    if (boost::iequals(property_type, "Buddenberg_Wilke_1949"))
    {
        DBUG("TODO: Implementation of %s property!!", property_type.c_str());
        return std::make_unique<Constant>(0);
    }
    if (boost::iequals(property_type, "Peng_Robinson_1976"))
    {
        DBUG("TODO: Implementation of %s property!!", property_type.c_str());
        return std::make_unique<Constant>(0);
    }
    if (boost::iequals(property_type, "Brooks_Corey_1964"))
    {
        DBUG("TODO: Implementation of %s property!!", property_type.c_str());
        return std::make_unique<Constant>(0);
    }
    // If none of the above property types are found, OGS throws an error.
    OGS_FATAL(
        "The specified component property type \"%s\" was not "
        "recognized",
        property_type.c_str());
    return nullptr;  // to avoid 'no return' warnings
}

}  // MaterialPropertyLib
