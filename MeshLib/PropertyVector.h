// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <cassert>
#include <iterator>
#include <range/v3/range/concepts.hpp>
#include <range/v3/view/common.hpp>
#include <string>
#ifdef _MSC_VER
#include <string_view>  // for MSVC, we need to include <string_view> to use std::string_view
#endif
#include <utility>
#include <vector>

#include "BaseLib/Algorithm.h"
#include "MeshEnums.h"

namespace MeshLib
{

#ifdef _MSC_VER  // MSVC does not support inline variables until C++17, so we
                 // use a constexpr string_view instead of an inline variable.
inline constexpr std::string_view vtkGhostTypeString = "vtkGhostType";
#else
inline constexpr std::string vtkGhostTypeString = "vtkGhostType";
#endif

class PropertyVectorBase
{
public:
    virtual PropertyVectorBase* clone(
        std::vector<std::size_t> const& exclude_positions) const = 0;
    virtual ~PropertyVectorBase() = default;

    MeshItemType getMeshItemType() const { return _mesh_item_type; }
    std::string const& getPropertyName() const { return _property_name; }
    int getNumberOfGlobalComponents() const { return _n_components; }

protected:
    PropertyVectorBase(std::string property_name,
                       MeshItemType mesh_item_type,
                       std::size_t n_components)
        : _n_components(n_components),
          _mesh_item_type(mesh_item_type),
          _property_name(std::move(property_name))
    {
    }

    int const _n_components;
    MeshItemType const _mesh_item_type;
    std::string const _property_name;
};

/// Contiguous storage for the values of one named mesh property.
/// Values are laid out tuple-contiguous: the component of a tuple is at
/// tuple_index * n_components + component, so the components of one tuple are
/// adjacent in memory. Instances are created and owned by Properties.
/// \tparam PROP_VAL_TYPE typically a scalar, a vector or a matrix
template <typename PROP_VAL_TYPE>
class PropertyVector : public PropertyVectorBase
{
    friend class Properties;

public:
    using value_type = PROP_VAL_TYPE;

public:
    constexpr std::size_t getNumberOfTuples() const
    {
        return size() / _n_components;
    }

    //! Returns the value for the given component stored in the given tuple.
    PROP_VAL_TYPE& getComponent(std::size_t tuple_index, int component)
    {
        assert(component < _n_components);
        assert(tuple_index < getNumberOfTuples());
        return data_[tuple_index * getNumberOfGlobalComponents() + component];
    }

    //! Returns the value for the given component stored in the given tuple.
    PROP_VAL_TYPE const& getComponent(std::size_t tuple_index,
                                      int component) const
    {
        assert(component < _n_components);
        assert(tuple_index < getNumberOfTuples());
        return data_[tuple_index * getNumberOfGlobalComponents() + component];
    }

    PropertyVectorBase* clone(
        std::vector<std::size_t> const& exclude_positions) const override
    {
        auto* cloned_pv = new PropertyVector<PROP_VAL_TYPE>(
            _property_name, _mesh_item_type, _n_components);
        cloned_pv->data_ = BaseLib::excludeObjectCopy(data_, exclude_positions);
        return cloned_pv;
    }

    /// Method returns the number of tuples times the number of tuple
    /// components.
    constexpr std::size_t size() const { return data_.size(); }

    constexpr std::ptrdiff_t ssize() const { return std::ssize(data_); }

    // Same as begin, but semantically different; begin, end are pairs of
    // iterators, data returns raw pointer.
    constexpr const PROP_VAL_TYPE* data() const { return data_.data(); }
    constexpr PROP_VAL_TYPE* data() { return data_.data(); }

    constexpr PROP_VAL_TYPE* begin() { return data_.data(); }
    constexpr PROP_VAL_TYPE* end() { return data_.data() + data_.size(); }

    constexpr const PROP_VAL_TYPE* cbegin() const { return data_.data(); }
    constexpr const PROP_VAL_TYPE* cend() const
    {
        return data_.data() + data_.size();
    }
    constexpr const PROP_VAL_TYPE* begin() const { return cbegin(); }
    constexpr const PROP_VAL_TYPE* end() const { return cend(); }

    constexpr PROP_VAL_TYPE& operator[](std::size_t const pos)
    {
        if constexpr (std::is_same_v<PROP_VAL_TYPE, bool>)
        {
            static_assert(!std::is_same_v<PROP_VAL_TYPE, bool>,
                          "PropertyVector<bool>::operator[] cannot be "
                          "instantiated for booleans.");
        }
        else
        {
            return data_[pos];
        }
    }
    constexpr PROP_VAL_TYPE const& operator[](std::size_t const pos) const
    {
        return data_[pos];
    }

    constexpr void resize(std::size_t const size) { data_.resize(size); }
    constexpr void resize(std::size_t const size, const PROP_VAL_TYPE& value)
    {
        data_.resize(size, value);
    }

    template <typename R>
        requires std::ranges::input_range<R> &&
                 std::convertible_to<std::ranges::range_value_t<R>,
                                     PROP_VAL_TYPE>
    constexpr void assign(R&& r)
    {
#if __cpp_lib_containers_ranges >= 202202L
        data_.assign_range(r);
#else
        if constexpr (ranges::common_range<R>)
        {
            data_.assign(r.begin(), r.end());
        }
        else
        {
            data_.assign(ranges::views::common(r).begin(),
                         ranges::views::common(r).end());
        }
#endif
    }

    constexpr void push_back(const PROP_VAL_TYPE& value)
    {
        data_.push_back(value);
    }

    constexpr void clear() { data_.clear(); }
    constexpr bool empty() const { return data_.empty(); }

protected:
    /// @brief The constructor taking meta information for the data.
    /// @param property_name a string describing the property
    /// @param mesh_item_type the values of the property are either assigned to
    /// nodes or cells (see enumeration MeshItemType)
    /// @param n_components the number of components of a property
    explicit PropertyVector(std::string const& property_name,
                            MeshItemType mesh_item_type,
                            std::size_t n_components)
        : PropertyVectorBase(property_name, mesh_item_type, n_components)
    {
    }

    /// @brief The constructor taking meta information for the data.
    /// @param n_property_values number of property values (value can be a tuple
    /// with several entries)
    /// @param property_name a string describing the property
    /// @param mesh_item_type the values of the property are either assigned to
    /// nodes or cells (see enumeration MeshItemType)
    /// @param n_components the number of components of a property
    PropertyVector(std::size_t n_property_values,
                   std::string const& property_name,
                   MeshItemType mesh_item_type,
                   std::size_t n_components)
        : PropertyVectorBase(property_name, mesh_item_type, n_components),
          data_(n_property_values * n_components)
    {
    }

private:
    std::vector<PROP_VAL_TYPE> data_;
};

static_assert(ranges::contiguous_range<PropertyVector<double>>);
static_assert(ranges::sized_range<PropertyVector<double>>);
}  // end namespace MeshLib
