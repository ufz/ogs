/**
 * \file
 * \author Norbert Grunwald
 * \date   Sep 7, 2017
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */
#pragma once

#include <array>
#include <boost/variant.hpp>
#include <string>

#include "BaseLib/ConfigTree.h"

#include "mpEnums.h"

namespace MaterialPropertyLib
{
class Medium;
class Phase;
class Component;

/// This is a custom data type for arbitrary properties, based on the
/// boost::variant container. It can hold scalars, vectors, tensors, and so
/// forth.
enum PropertyDataTypeName
{
    nScalar,
    nPair,
    nVector,
    nSymmTensor,
    nTensor
};

using PropertyDataType =
    boost::variant<double, Pair, Vector, SymmTensor, Tensor, std::string>;

/// This class is the base class for any material property of any
/// scale (i.e. components, phases, media, ...). The single value of
/// that Property can hold scalars, vectors, tensors, strings, etc.
class Property
{
protected:
    /// The single value of a property.
    PropertyDataType _value;
    PropertyDataType _dvalue;

public:
    virtual ~Property() = default;
    /// This method is called when a property is used for the wrong kind of
    /// material, or if the property is not implemented on this kind of material
    /// yet.
    void notImplemented(std::string /*property*/, std::string /*material*/);
    /// This virtual method simply returns the private _value attribute without
    /// changing it.
    virtual PropertyDataType value() const;
    /// This virtual method will compute the property value based on the primary
    /// variables that are passed as arguments.
    virtual PropertyDataType value(VariableArray const& /*unused*/);
    /// This virtual method will compute the derivative of a property
    /// with respect to the given variables pv.
    virtual PropertyDataType dvalue(VariableArray const&, Variables const pv);
    /// This virtual method will compute the second derivative of a
    /// property with respect to the given variables pv1 and pv2.
    virtual PropertyDataType ddvalue(VariableArray const&, Variables const pv1,
                                     Variables const pv2);
};

/// This data type is based on a std::array. It can hold pointers to objects of
/// class Property or its inheritors. The size of this array is determined by
/// the number of entries of the PropertyEnum enumerator.
using PropertyArray =
    std::array<std::unique_ptr<Property>, number_of_property_enums>;

/// Method to select a property by name and to call a derived property
/// constructor.
///
/// This method creates a new Property. It uses the information stored inside
/// the configuration tree to select the property that is specified by the user.
/// It further passes a pointer to the material that requests the property.
template <typename MaterialType>
std::unique_ptr<Property> selectProperty(BaseLib::ConfigTree const& config,
                                         MaterialType material_type);
/// This method creates a new medium property.
std::unique_ptr<Property> newProperty(BaseLib::ConfigTree const& config,
                                      Medium* m);
/// This method creates a new phase property.
std::unique_ptr<Property> newProperty(BaseLib::ConfigTree const& config,
                                      Phase* p);
/// This method creates a new component property.
std::unique_ptr<Property> newProperty(BaseLib::ConfigTree const& config,
                                      Component* c);
/// This method returns the 0-based index of the variant data types. Can be
/// enhanced by using enums.
inline std::size_t getType(Property const& p)
{
    return p.value().which();
}

/// This method returns a value of type double from the property value
/// attribute.
inline double getScalar(Property& p)
{
    assert((getType(p) == PropertyDataTypeName::nScalar) &&
           "The requested "
           "value type is not of type 'double'");

    return boost::get<double>(p.value());
}

/// This method forces the computation of a value of type double and returns it.
inline double getScalar(Property& p, VariableArray const& v)
{
    return boost::get<double>(p.value(v));
}

/// This method forces the computation of a value of type Pair and returns it
inline Pair getPair(Property& p, VariableArray const& v)
{
    return boost::get<Pair>(p.value(v));
}

/// This method forces the computation of a value of type Vector and returns it.
inline Vector getVector(Property& p)
{
    assert((getType(p) == PropertyDataTypeName::nVector) &&
           "The requested value type is not of type 'Vector'");
    return boost::get<Vector>(p.value());
}

/// This method forces the computation of a value of type Vector and returns it.
inline Vector getVector(Property& p, VariableArray const& v)
{
    return boost::get<Vector>(p.value(v));
}

/// This method forces the computation the derivative of a value, returns a pair
/// of double values. The derivative is computed with respect to variable pv.
inline Pair getPairDerivative(Property& p, VariableArray const& v,
                              Variables const pv)
{
    return boost::get<Pair>(p.dvalue(v, pv));
}

/// This method forces the computation the derivative of a value, returns a
/// double value. The derivative is computed with respect to variable pv.
inline double getScalarDerivative(Property& p, VariableArray const& v,
                                  Variables const pv)
{
    return boost::get<double>(p.dvalue(v, pv));
}

/// This method forces the computation the second derivative of a value, returns
/// a double value. The derivative is computed with respect to variables pv1 and
/// pv2.
inline double getScalarDerivative(Property& p, VariableArray const& v,
                                  Variables const pv1, Variables const pv2)
{
    return boost::get<double>(p.ddvalue(v, pv1, pv2));
}

/// This method returns a value of type string from the property value
/// attribute.
inline std::string getString(Property const& p)
{
    return boost::get<std::string>(p.value());
}

/// This method returns a value of any valid type from the property value
/// attribute. The data type is provided by the first parameter in the argument
/// list.
template <typename T>
T getValue(T const&, Property const& p)
{
    return boost::get<T>(p.value());
}

}  // namespace MaterialPropertyLib
