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
#pragma once

#include <array>
#include <string>
#include <variant>

#include "PropertyType.h"
#include "VariableType.h"

#include "ParameterLib/SpatialPosition.h"

namespace MaterialPropertyLib
{
/// This is a custom data type for arbitrary properties, based on the
/// std::variant container. It can hold scalars, vectors, tensors, and so
/// forth.
enum PropertyDataTypeName
{
    nScalar,
    nPair,
    nVector,
    nSymmTensor,
    nTensor
};

using PropertyDataType = std::
    variant<double, Pair, Vector, Tensor2d, SymmTensor, Tensor, std::string>;

/// This class is the base class for any material property of any
/// scale (i.e. components, phases, media, ...). The single value of
/// that Property can hold scalars, vectors, tensors, strings, etc.
class Property
{
public:
    virtual ~Property() = default;
    /// This method is called when a property is used for the wrong kind of
    /// material, or if the property is not implemented on this kind of material
    /// yet.
    void notImplemented(const std::string& property,
                        const std::string& material) const;
    /// This virtual method simply returns the private _value attribute without
    /// changing it.
    virtual PropertyDataType value() const;
    /// This virtual method will compute the property value based on the primary
    /// variables that are passed as arguments.
    virtual PropertyDataType value(VariableArray const& variable_array,
                                   ParameterLib::SpatialPosition const& pos,
                                   double const t) const;
    /// This virtual method will compute the derivative of a property
    /// with respect to the given variable pv.
    virtual PropertyDataType dValue(VariableArray const& variable_array,
                                    Variable const variable) const;
    /// This virtual method will compute the second derivative of a
    /// property with respect to the given variables pv1 and pv2.
    virtual PropertyDataType d2Value(VariableArray const& variable_array,
                                     Variable const variable1,
                                     Variable const variable2) const;

    template <typename T>
    T value() const
    {
        return std::get<T>(value());
    }

    template <typename T>
    T value(VariableArray const& variable_array,
            ParameterLib::SpatialPosition const& pos,
            double const t) const
    {
        return std::get<T>(value(variable_array, pos, t));
    }

    template <typename T>
    T dValue(VariableArray const& variable_array,
             Variable const variable) const
    {
        return std::get<T>(dValue(variable_array, variable));
    }
    template <typename T>
    T d2Value(VariableArray const& variable_array,
              Variable const& variable1,
              Variable const& variable2) const
    {
        return std::get<T>(d2Value(variable_array, variable1, variable2));
    }

protected:
    /// The single value of a property.
    PropertyDataType _value;
    PropertyDataType _dvalue;
};

/// This method returns the 0-based index of the variant data types. Can be
/// enhanced by using enums.
inline std::size_t getType(Property const& p)
{
    return p.value().index();
}

inline void overwriteExistingProperties(PropertyArray& properties,
                                        PropertyArray& new_properties)
{
    for (std::size_t i = 0; i < properties.size(); ++i)
    {
        if (new_properties[i] != nullptr)
        {
            properties[i] = std::move(new_properties[i]);
        }
    }
}

}  // namespace MaterialPropertyLib
