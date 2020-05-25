/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <cassert>
#include <iterator>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

#include "BaseLib/Algorithm.h"
#include "Location.h"

namespace MeshLib
{
class PropertyVectorBase
{
public:
    virtual PropertyVectorBase* clone(
        std::vector<std::size_t> const& exclude_positions) const = 0;
    virtual ~PropertyVectorBase() = default;

    MeshItemType getMeshItemType() const { return mesh_item_type_; }
    std::string const& getPropertyName() const { return property_name_; }
    int getNumberOfComponents() const { return n_components_; }

protected:
    PropertyVectorBase(std::string property_name,
                       MeshItemType mesh_item_type,
                       std::size_t n_components)
        : n_components_(n_components),
          mesh_item_type_(mesh_item_type),
          property_name_(std::move(property_name))
    {
    }

    int const n_components_;
    MeshItemType const mesh_item_type_;
    std::string const property_name_;
};

/// Class template PropertyVector is a std::vector with template parameter
/// PROP_VAL_TYPE. The reason for the derivation of std::vector is
/// the template specialisation for pointer types below.
/// \tparam PROP_VAL_TYPE typical this is a scalar, a vector or a matrix
template <typename PROP_VAL_TYPE>
class PropertyVector : public std::vector<PROP_VAL_TYPE>,
                       public PropertyVectorBase
{
    friend class Properties;

public:
    std::size_t getNumberOfTuples() const
    {
        return std::vector<PROP_VAL_TYPE>::size() / n_components_;
    }

    //! Returns the value for the given component stored in the given tuple.
    PROP_VAL_TYPE& getComponent(std::size_t tuple_index, int component)
    {
        assert(component < n_components_);
        assert(tuple_index < getNumberOfTuples());
        return this->operator[](tuple_index* getNumberOfComponents() +
                                component);
    }

    //! Returns the value for the given component stored in the given tuple.
    PROP_VAL_TYPE const& getComponent(std::size_t tuple_index,
                                      int component) const
    {
        assert(component < n_components_);
        assert(tuple_index < getNumberOfTuples());
        return this->operator[](tuple_index* getNumberOfComponents() +
                                component);
    }

    PropertyVectorBase* clone(
        std::vector<std::size_t> const& exclude_positions) const override
    {
        auto* t(new PropertyVector<PROP_VAL_TYPE>(
            property_name_, mesh_item_type_, n_components_));
        BaseLib::excludeObjectCopy(*this, exclude_positions, *t);
        return t;
    }

    /// Method returns the number of tuples times the number of tuple
    /// components.
    std::size_t size() const { return std::vector<PROP_VAL_TYPE>::size(); }

protected:
    /// @brief The constructor taking meta information for the data.
    /// @param property_name a string describing the property
    /// @param mesh_item_type the values of the property are either assigned to
    /// nodes or cells (see enumeration MeshItemType)
    /// @param n_components the number of components of a property
    explicit PropertyVector(std::string const& property_name,
                            MeshItemType mesh_item_type,
                            std::size_t n_components)
        : std::vector<PROP_VAL_TYPE>(),
          PropertyVectorBase(property_name, mesh_item_type, n_components)
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
        : std::vector<PROP_VAL_TYPE>(n_property_values * n_components),
          PropertyVectorBase(property_name, mesh_item_type, n_components)
    {
    }
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
class PropertyVector<T*> : public std::vector<std::size_t>,
                           public PropertyVectorBase
{
    friend class Properties;

public:
    /// Destructor ensures the deletion of the heap-constructed objects.
    ~PropertyVector() override
    {
        for (auto v : values_)
        {
            delete[] v;
        }
    }

    /// The operator[] uses the item to group property map to access to the
    /// correct property value/object.
    T* const& operator[](std::size_t id) const
    {
        return values_[std::vector<std::size_t>::operator[](id)];
    }

    T*& operator[](std::size_t id)
    {
        return values_[std::vector<std::size_t>::operator[](id)];
    }

    void initPropertyValue(std::size_t group_id, T const& value)
    {
        if (n_components_ != 1)
        {
            OGS_FATAL(
                "Single-component version of initPropertyValue() is called "
                "for a multi-components PropertyVector<T*>");
        }
        auto* p = new T[1];
        p[0] = value;
        values_[group_id] = p;
    }

    void initPropertyValue(std::size_t group_id, std::vector<T> const& values)
    {
        if (n_components_ != static_cast<int>(values.size()))
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
        values_[group_id] = p;
    }

    std::size_t getNumberOfTuples() const
    {
        return std::vector<std::size_t>::size();
    }

    /// Method returns the number of tuples times the number of tuple
    /// components.
    std::size_t size() const
    {
        return n_components_ * std::vector<std::size_t>::size();
    }

    PropertyVectorBase* clone(
        std::vector<std::size_t> const& exclude_positions) const override
    {
        // create new PropertyVector with modified mapping
        PropertyVector<T*>* t(new PropertyVector<T*>(
            values_.size() / n_components_,
            BaseLib::excludeObjectCopy(*this, exclude_positions),
            property_name_, mesh_item_type_, n_components_));
        // copy pointers to property values
        for (std::size_t j(0); j < values_.size(); j++)
        {
            std::vector<T> values(values_[j], values_[j] + n_components_);
            t->initPropertyValue(j, values);
        }
        return t;
    }

    //! Returns the value for the given component stored in the given tuple.
    T const& getComponent(std::size_t tuple_index, int component) const
    {
        assert(component < n_components_);
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

#ifndef NDEBUG
    std::ostream& print(std::ostream& os) const
    {
        os << "\nPropertyVector<T*> at address: " << this << ":\n";
        os << "\tmapping (" << size() << "):\n";
        std::copy(this->cbegin(), this->cend(),
                  std::ostream_iterator<std::size_t>(os, " "));
        os << "\n\tvalues (" << values_.size() << "):\n";
        for (std::size_t k(0); k < values_.size(); k++)
        {
            os << "val: " << *(values_[k]) << ", address: " << values_[k]
               << "\n";
        }
        return os;
    }
#endif

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
        : std::vector<std::size_t>(std::move(item2group_mapping)),
          PropertyVectorBase(property_name, mesh_item_type, n_components),
          values_(n_prop_groups * n_components)
    {
    }

private:
    std::vector<T*> values_;
    // hide method
    T* at(std::size_t);
};

}  // end namespace MeshLib
