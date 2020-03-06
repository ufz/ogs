/**
 * \file
 * \author Norbert Grunwald
 * \date   Sep 7, 2017
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */
#pragma once

#include <Eigen/Dense>
#include <array>
#include <string>
#include <variant>

#include "ParameterLib/SpatialPosition.h"
#include "PropertyType.h"
#include "VariableType.h"

namespace MaterialPropertyLib
{
class Medium;
class Phase;
class Component;

using PropertyDataType =
    std::variant<double, Eigen::Matrix<double, 2, 1>,
                 Eigen::Matrix<double, 3, 1>, Eigen::Matrix<double, 2, 2>,
                 Eigen::Matrix<double, 3, 3>, Eigen::Matrix<double, 4, 1>,
                 Eigen::Matrix<double, 6, 1>>;

/// Conversion of a vector to PropertyDataType for different sizes of the
/// vector.
/// \attention It cannot distinguish between 2x2 matrix and 4x1 vector.
PropertyDataType fromVector(std::vector<double> const& values);

/// This class is the base class for any material property of any
/// scale (i.e. components, phases, media, ...). The single value of
/// that Property can hold scalars, vectors, tensors, strings, etc.
class Property
{
public:
    virtual ~Property() = default;

    /// Returns the initial (or reference) value of the property.
    /// The default implementation forwards to the value function.
    virtual PropertyDataType initialValue(
        ParameterLib::SpatialPosition const& pos, double const t) const;

    /// This virtual method simply returns the private _value attribute without
    /// changing it.
    virtual PropertyDataType value() const;
    /// This virtual method will compute the property value based on the primary
    /// variables that are passed as arguments.
    virtual PropertyDataType value(VariableArray const& variable_array,
                                   ParameterLib::SpatialPosition const& pos,
                                   double const t, double const dt) const;
    /// This virtual method will compute the derivative of a property
    /// with respect to the given variable pv.
    virtual PropertyDataType dValue(VariableArray const& variable_array,
                                    Variable const variable,
                                    ParameterLib::SpatialPosition const& pos,
                                    double const t, double const dt) const;
    /// This virtual method will compute the second derivative of a
    /// property with respect to the given variables pv1 and pv2.
    virtual PropertyDataType d2Value(VariableArray const& variable_array,
                                     Variable const variable1,
                                     Variable const variable2,
                                     ParameterLib::SpatialPosition const& pos,
                                     double const t, double const dt) const;

    /// This virtual method will compute the inverse value of the property
    /// based on the primary variables that are passed as arguments.
    virtual PropertyDataType inverse_value(
        VariableArray const& variable_array,
        ParameterLib::SpatialPosition const& pos, double const t,
        double const dt) const;

    virtual void setScale(
        std::variant<Medium*, Phase*, Component*> /*scale_pointer*/){};

    template <typename T>
    T initialValue(ParameterLib::SpatialPosition const& pos,
                   double const t) const
    {
        return std::get<T>(initialValue(pos, t));
    }

    template <typename T>
    T value() const
    {
        return std::get<T>(value());
    }

    template <typename T>
    T value(VariableArray const& variable_array,
            ParameterLib::SpatialPosition const& pos, double const t,
            double const dt) const
    {
        return std::get<T>(value(variable_array, pos, t, dt));
    }

    template <typename T>
    T inverse_value(VariableArray const& variable_array,
                    ParameterLib::SpatialPosition const& pos, double const t,
                    double const dt) const
    {
        return std::get<T>(inverse_value(variable_array, pos, t, dt));
    }

    template <typename T>
    T dValue(VariableArray const& variable_array, Variable const variable,
             ParameterLib::SpatialPosition const& pos, double const t,
             double const dt) const
    {
        return std::get<T>(dValue(variable_array, variable, pos, t, dt));
    }
    template <typename T>
    T d2Value(VariableArray const& variable_array, Variable const& variable1,
              Variable const& variable2,
              ParameterLib::SpatialPosition const& pos, double const t,
              double const dt) const
    {
        return std::get<T>(
            d2Value(variable_array, variable1, variable2, pos, t, dt));
    }

protected:
    /// The single value of a property.
    PropertyDataType _value;
    PropertyDataType _dvalue;
};

inline void overwriteExistingProperties(
    PropertyArray& properties,
    PropertyArray& new_properties,
    std::variant<Medium*, Phase*, Component*>
        scale_pointer)
{
    for (std::size_t i = 0; i < properties.size(); ++i)
    {
        if (new_properties[i] != nullptr)
        {
            properties[i] = std::move(new_properties[i]);
            properties[i]->setScale(scale_pointer);
        }
    }
}

}  // namespace MaterialPropertyLib
