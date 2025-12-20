// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <cassert>
#include <iterator>
#include <range/v3/range/concepts.hpp>
#include <range/v3/view/common.hpp>
#include <string>
#include <utility>
#include <vector>

#include "BaseLib/Algorithm.h"
#include "MeshEnums.h"

namespace MeshLib
{
class PropertyVectorBase
{
public:
    virtual PropertyVectorBase* clone(
        std::vector<std::size_t> const& exclude_positions) const = 0;
    virtual ~PropertyVectorBase() = default;

    MeshItemType getMeshItemType() const { return _mesh_item_type; }
    std::string const& getPropertyName() const { return _property_name; }
    int getNumberOfGlobalComponents() const { return _n_components; }
    bool is_for_output = true;

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

/// Class template PropertyVector is a std::vector with template parameter
/// PROP_VAL_TYPE. The reason for the derivation of std::vector is
/// the template specialisation for pointer types below.
/// \tparam PROP_VAL_TYPE typical this is a scalar, a vector or a matrix
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

/// Class template PropertyVector is a std::vector with template parameter
/// T, where T is a pointer type.
/// The behaviour has changed for the constructor, destructor and the
/// operator[]. The user has to provide the size and an item to group mapping
/// for construction. The destructor takes care to delete the entries of the
/// vector. The operator[] uses an item-to-group property map to access the
/// correct property.
/// \tparam T pointer type, the type the type points to is typical a scalar,
/// a vector or a matrix type
template <typename T>
class PropertyVector<T*> : public PropertyVectorBase
{
    friend class Properties;

public:
    /// Destructor ensures the deletion of the heap-constructed objects.
    ~PropertyVector() override
    {
        for (auto v : _values)
        {
            delete[] v;
        }
    }

    /// The operator[] uses the item to group property map to access to the
    /// correct property value/object.
    T* const& operator[](std::size_t id) const
    {
        return _values[_item2group_mapping[id]];
    }

    T*& operator[](std::size_t id) { return _values[_item2group_mapping[id]]; }

    void initPropertyValue(std::size_t group_id, T const& value)
    {
        if (_n_components != 1)
        {
            OGS_FATAL(
                "Single-component version of initPropertyValue() is called "
                "for a multi-components PropertyVector<T*>");
        }
        auto* p = new T[1];
        p[0] = value;
        _values[group_id] = p;
    }

    void initPropertyValue(std::size_t group_id, std::vector<T> const& values)
    {
        if (_n_components != static_cast<int>(values.size()))
        {
            OGS_FATAL(
                "The size of provided values in initPropertyValue() is "
                "not same as the number of components in PropertyVector<T*>");
        }

        auto* p = new T[values.size()];
        for (unsigned i = 0; i < values.size(); i++)
        {
            p[i] = values[i];
        }
        _values[group_id] = p;
    }

    std::size_t getNumberOfTuples() const { return _item2group_mapping.size(); }

    /// Method returns the number of tuples times the number of tuple
    /// components.
    std::size_t size() const { return _n_components * getNumberOfTuples(); }

    PropertyVectorBase* clone(
        std::vector<std::size_t> const& exclude_positions) const override
    {
        // create new PropertyVector with modified mapping
        PropertyVector<T*>* t(new PropertyVector<T*>(
            _values.size() / _n_components,
            BaseLib::excludeObjectCopy(_item2group_mapping, exclude_positions),
            _property_name, _mesh_item_type, _n_components));
        // copy pointers to property values
        for (std::size_t j(0); j < _values.size(); j++)
        {
            std::vector<T> values(_values[j], _values[j] + _n_components);
            t->initPropertyValue(j, values);
        }
        return t;
    }

    //! Returns the value for the given component stored in the given tuple.
    T const& getComponent(std::size_t tuple_index, int component) const
    {
        assert(component < _n_components);
        assert(tuple_index < getNumberOfTuples());

        const double* p = this->operator[](tuple_index);
        if (p == nullptr)
        {
            OGS_FATAL(
                "No data found in the property vector {:s} "
                "for the tuple index {:d} and component {:d}",
                getPropertyName(), tuple_index, component);
        }
        return p[component];
    }

protected:
    /// @brief The constructor taking meta information for the data.
    /// @param n_prop_groups number of different property values
    /// @param item2group_mapping Class Mesh has a mapping from the mesh items
    /// (Node or Element) to an index (position in the data structure).
    /// The vector item2group_mapping must have the same number of entries as
    /// the above mapping and the values have to be in the range
    /// \f$[0, \text{n\_prop\_groups})\f$.
    /// @param property_name a string describing the property
    /// @param mesh_item_type the values of the property are either assigned to
    /// nodes or cells (see enumeration MeshItemType)
    /// @param n_components the number of elements of a tuple
    PropertyVector(std::size_t n_prop_groups,
                   std::vector<std::size_t>
                       item2group_mapping,
                   std::string const& property_name,
                   MeshItemType mesh_item_type,
                   std::size_t n_components)
        : PropertyVectorBase(property_name, mesh_item_type, n_components),
          _item2group_mapping(std::move(item2group_mapping)),
          _values(n_prop_groups * n_components)
    {
    }

private:
    std::vector<std::size_t> _item2group_mapping;
    std::vector<T*> _values;
    // hide method
    T* at(std::size_t);
};

static_assert(ranges::contiguous_range<PropertyVector<double>>);
static_assert(ranges::sized_range<PropertyVector<double>>);
}  // end namespace MeshLib
