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
#ifndef MATERIALLIB_MPL_MPPROPERTY_H_
#define MATERIALLIB_MPL_MPPROPERTY_H_

#include "BaseLib/ConfigTree.h"
#include "mpEnums.h"
#include <array>
#include <boost/variant.hpp>
#include <string>

namespace MaterialPropertyLib
{

class Medium;
class Phase;
class Component;

/**
 * This is a custom data type for arbitrary properties, based on the
 * boost::variant container. It can hold scalars, vectors, tensors, and
 * strings. It can be further extended to hold symmetrical tensor data,
 * which consist of six components.
*/
using PropertyDataType = boost::variant<double, Vector, Tensor, std::string>;

/**
 * This class is the base class for any material property of any
 * scale (i.e. components, phases, media, ...). The single value of
 * that Property can hold scalars, vectors, tensors, or strings
 */
class Property
{
protected:
	/// The single value of a property.
    PropertyDataType _value;
    bool _set;
public:
    Property();
    void set(bool);
    bool set ();
    virtual ~Property(){};
    /// This method is called when a property is used for the wrong
    /// kind of material, or if the proeprty is not implemented on
    /// this kind of material yet.
    void notImplemented (std::string, std::string);
    /// This virtual method simply returns the private _value attribute
    /// without changing it.
    virtual PropertyDataType value () const;
    /// This virtual method will compute the property value baesd on the
    /// primary variables that are passed as arguments.
    virtual PropertyDataType value (VariableArray const&);
};  // class Property

/**
 * This data type is based on a std::array. It can hold pointers
 * to objects of class Property or its inheritors. The size of
 * this array is determined by the number of entries of the
 * PropertyEnum enumerator.
*/
using PropertyArray = std::array<std::unique_ptr<Property>, number_of_property_enums>;

/// Method to select a property by name and to call a derived property
/// constructor.
template <typename MaterialType>
//Property* selectProperty (BaseLib::ConfigTree const&, MaterialType);
std::unique_ptr<Property> selectProperty
(BaseLib::ConfigTree const&, MaterialType);
/// This method creates a new medium property.
std::unique_ptr<Property> newProperty(BaseLib::ConfigTree const& config, Medium*);
/// This method creates a new phase property.
std::unique_ptr<Property> newProperty(BaseLib::ConfigTree const&, Phase*);
/// This method creates a new component property.
std::unique_ptr<Property> newProperty(BaseLib::ConfigTree const& config, Component*);

/// This method returns a value of type double from the
/// property value attribute
inline double getScalar (Property* p)
{
    return boost::get<double>(p->value());
}

/// This method forces the computation of a value of type double
/// and returns it
inline double getScalar (Property* p, VariableArray const&v)
{
    return boost::get<double>(p->value(v));
}

/// This method returns a value of type string from the
/// property value attribute
inline std::string getString (Property* p)
{
    return boost::get<std::string>(p->value());
}
/// This method returns a value of any valid type from
/// the property value attribute. The data type is provided
/// by the first parameter in the argument list.
template <typename T>
T getValue (T const&, Property* p)
{
    return boost::get<T>(p->value());
}

} //MaterialPropertyLib


#endif /* MATERIALLIB_MPL_MPPROPERTY_H_ */
