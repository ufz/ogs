/**
 * \file
 * \author Norbert Grunwald
 * \date   Sep 7, 2017
 *
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */
#pragma once

#include <Eigen/Dense>
#include <array>
#include <string>
#include <typeinfo>
#include <variant>

#include "BaseLib/Error.h"
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
                 Eigen::Matrix<double, 6, 1>, Eigen::MatrixXd>;

/// Conversion of a vector to PropertyDataType for different sizes of the
/// vector.
///
/// \attention This method cannot distinguish between 2x2 matrix and 4x1 vector.
///
/// \note If the \c values vector stores all elements of a 2x2 or 3x3 matrix
/// (i.e. the general case for 2x2 and 3x3 matrices), it is assumed to be in row
/// major storage order.
PropertyDataType fromVector(std::vector<double> const& values);

/// This class is the base class for any material property of any
/// scale (i.e. components, phases, media, ...). The single value of
/// that Property can hold scalars, vectors, tensors, strings, etc.
class Property
{
public:
#ifndef NDEBUG
    virtual ~Property()
    {
        if(property_used)
        {
            DBUG("Property is used: '{:s}'", description());
        }
        else
        {
            WARN("Property is not used: '{:s}'", description());
        }
    }
#else
    virtual ~Property() = default;
#endif

    /// Returns the initial (or reference) value of the property.
    /// The default implementation forwards to the value function.
    virtual PropertyDataType initialValue(
        ParameterLib::SpatialPosition const& pos, double const t) const;

    /// This virtual method simply returns the private value_ attribute without
    /// changing it.
    virtual PropertyDataType value() const;
    /// This virtual method will compute the property value based on the
    /// variables that are passed as arguments and the variables from the
    /// previous time step.
    virtual PropertyDataType value(VariableArray const& variable_array,
                                   VariableArray const& variable_array_prev,
                                   ParameterLib::SpatialPosition const& pos,
                                   double const t, double const dt) const;
    /// This virtual method will compute the property value based on the
    /// variables that are passed as arguments with the default implementation
    /// using empty variables array for the previous time step.
    virtual PropertyDataType value(VariableArray const& variable_array,
                                   ParameterLib::SpatialPosition const& pos,
                                   double const t, double const dt) const;
    /// This virtual method will compute the property derivative value based on
    /// the variables that are passed as arguments and the variables from the
    /// previous time step.
    virtual PropertyDataType dValue(VariableArray const& variable_array,
                                    VariableArray const& variable_array_prev,
                                    Variable const variable,
                                    ParameterLib::SpatialPosition const& pos,
                                    double const t, double const dt) const;
    /// This virtual method will compute the property derivative value based on
    /// the variables that are passed as arguments with the default
    /// implementation using empty variables array for the previous time step.
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

    void setScale(std::variant<Medium*, Phase*, Component*> scale)
    {
        scale_ = scale;
        checkScale();
    };

    template <typename T>
    T initialValue(ParameterLib::SpatialPosition const& pos,
                   double const t) const
    {
        try
        {
            return std::get<T>(initialValue(pos, t));
        }
        catch (std::bad_variant_access const&)
        {
            OGS_FATAL(
                "The initial value of {:s} does not hold requested type '{:s}' "
                "but a {:s}.",
                description(),
                typeid(T).name(),
                property_data_type_names_[initialValue(pos, t).index()]);
        }
    }

    template <typename T>
    T value() const
    {
        try
        {
#ifndef NDEBUG
            property_used = true;
#endif
            return std::get<T>(value());
        }
        catch (std::bad_variant_access const&)
        {
            OGS_FATAL(
                "The value of {:s} does not hold requested type '{:s}' but a "
                "{:s}.",
                description(),
                typeid(T).name(),
                property_data_type_names_[value().index()]);
        }
    }

    template <typename T>
    T value(VariableArray const& variable_array,
            VariableArray const& variable_array_prev,
            ParameterLib::SpatialPosition const& pos, double const t,
            double const dt) const
    {
        try
        {
#ifndef NDEBUG
            property_used = true;
#endif
            return std::get<T>(
                value(variable_array, variable_array_prev, pos, t, dt));
        }
        catch (std::bad_variant_access const&)
        {
            OGS_FATAL(
                "The value of {:s} is not of the requested type '{:s}' but a "
                "{:s}.",
                description(),
                typeid(T).name(),
                property_data_type_names_[value(variable_array,
                                                variable_array_prev, pos, t, dt)
                                              .index()]);
        }
    }
    template <typename T>
    T value(VariableArray const& variable_array,
            ParameterLib::SpatialPosition const& pos, double const t,
            double const dt) const
    {
        try
        {
#ifndef NDEBUG
            property_used = true;
#endif
            return std::get<T>(value(variable_array, pos, t, dt));
        }
        catch (std::bad_variant_access const&)
        {
            OGS_FATAL(
                "The value of {:s} is not of the requested type '{:s}' but a "
                "{:s}.",
                description(),
                typeid(T).name(),
                property_data_type_names_[value(variable_array, pos, t, dt)
                                              .index()]);
        }
    }

    template <typename T>
    T dValue(VariableArray const& variable_array,
             VariableArray const& variable_array_prev, Variable const variable,
             ParameterLib::SpatialPosition const& pos, double const t,
             double const dt) const
    {
        try
        {
#ifndef NDEBUG
            property_used = true;
#endif
            return std::get<T>(dValue(variable_array, variable_array_prev,
                                      variable, pos, t, dt));
        }
        catch (std::bad_variant_access const&)
        {
            OGS_FATAL(
                "The first derivative value of {:s} is not of the requested "
                "type '{:s}' but a {:s}.",
                description(),
                typeid(T).name(),
                property_data_type_names_
                    [dValue(variable_array, variable, pos, t, dt).index()]);
        }
    }
    template <typename T>
    T dValue(VariableArray const& variable_array, Variable const variable,
             ParameterLib::SpatialPosition const& pos, double const t,
             double const dt) const
    {
        try
        {
#ifndef NDEBUG
            property_used = true;
#endif
            return std::get<T>(dValue(variable_array, variable, pos, t, dt));
        }
        catch (std::bad_variant_access const&)
        {
            OGS_FATAL(
                "The first derivative value of {:s} is not of the requested "
                "type '{:s}' but a {:s}.",
                description(),
                typeid(T).name(),
                property_data_type_names_
                    [dValue(variable_array, variable, pos, t, dt).index()]);
        }
    }
    template <typename T>
    T d2Value(VariableArray const& variable_array, Variable const& variable1,
              Variable const& variable2,
              ParameterLib::SpatialPosition const& pos, double const t,
              double const dt) const
    {
        try
        {
#ifndef NDEBUG
            property_used = true;
#endif
            return std::get<T>(
                d2Value(variable_array, variable1, variable2, pos, t, dt));
        }
        catch (std::bad_variant_access const&)
        {
            OGS_FATAL(
                "The second derivative value of {:s} is not of the requested "
                "type '{:s}' but a {:s}.",
                description(),
                typeid(T).name(),
                property_data_type_names_[d2Value(variable_array, variable1,
                                                  variable2, pos, t, dt)
                                              .index()]);
        }
    }

protected:
    std::string name_;
    /// The single value of a property.
    PropertyDataType value_;
    PropertyDataType dvalue_;
    /// Definition scale of the property. Can be one of medium, phase, or
    /// component in general. Set through setScale method which takes care of
    /// the correctness in special cases.
    std::variant<Medium*, Phase*, Component*> scale_;

private:
    virtual void checkScale() const
    {
        // Empty check for properties which can be defined on every scale,
        // medium, phase or component
    }
    std::string description() const;
#ifndef NDEBUG
    mutable bool property_used = false;
#endif

private:
    /// Corresponds to the PropertyDataType
    static constexpr std::array property_data_type_names_ = {
        "scalar",           "2-vector",           "3-vector",
        "2x2-matrix",       "3x3-matrix",         "2D-Kelvin vector",
        "3D-Kelvin vector", "dynamic matrix type"};
    static_assert(property_data_type_names_.size() ==
                      std::variant_size_v<PropertyDataType>,
                  "The array of property data type names has different size "
                  "than the PropertyDataType variant type.");
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
