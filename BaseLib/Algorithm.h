/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <algorithm>
#include <cassert>
#include <concepts>
#include <optional>
#include <range/v3/algorithm/find_if.hpp>
#include <range/v3/range/concepts.hpp>
#include <range/v3/range/conversion.hpp>
#include <range/v3/view/concat.hpp>
#include <range/v3/view/partial_sum.hpp>
#include <range/v3/view/single.hpp>
#include <string>
#include <typeindex>
#include <typeinfo>
#include <utility>

#include "CompilerWorkarounds.h"
#include "Error.h"

namespace BaseLib
{
/// excludeObjectCopy copies only those objects that position within the source
/// vector is not in the exclude_positions vector. The implementation of the
/// algorithm requires that the given positions in exclude_positions are sorted
/// in ascending order.
/// @param src_vec the vector of source objects
/// @param exclude_positions the positions of objects in the source vector that
/// do not have to be copied
/// @return vector that contains the copied objects
template <typename T>
std::vector<T> excludeObjectCopy(
    std::vector<T> const& src_vec,
    std::vector<std::size_t> const& exclude_positions)
{
    std::vector<T> dest_vec;
    if (exclude_positions.empty())
    {
        dest_vec = src_vec;
        return dest_vec;
    }

    assert(exclude_positions.back() < src_vec.size());

    std::copy_n(src_vec.cbegin(), exclude_positions[0],
                std::back_inserter(dest_vec));
    for (std::size_t i = 1; i < exclude_positions.size(); ++i)
    {
        std::copy_n(src_vec.cbegin() + exclude_positions[i - 1] + 1,
                    exclude_positions[i] - (exclude_positions[i - 1] + 1),
                    std::back_inserter(dest_vec));
    }
    std::copy(src_vec.cbegin() + exclude_positions.back() + 1, src_vec.cend(),
              std::back_inserter(dest_vec));

    return dest_vec;
}

template <typename T>
void excludeObjectCopy(std::vector<T> const& src_vec,
                       std::vector<std::size_t> const& exclude_positions,
                       std::vector<T>& dest_vec)
{
    dest_vec = excludeObjectCopy(src_vec, exclude_positions);
}

/// Returns reference to an element in the range satisfying the predicate. If no
/// such element is found, error_callback is called and reference to
/// past-the-end of the range is returned.
template <ranges::input_range Range>
ranges::range_reference_t<Range> findElementOrError(
    Range& range,
    std::predicate<ranges::range_reference_t<Range>> auto&& predicate,
    std::invocable auto error_callback)
{
    auto it =
        ranges::find_if(range, std::forward<decltype(predicate)>(predicate));
    if (it == ranges::end(range))
    {
        error_callback();
        OGS_FATAL(
            "Element not found in the input range. The user provided error "
            "callback is meant not to return. That has not happened.");
    }
    return *it;
}

//! Inserts the given \c key with the given \c value into the \c map if an entry
//! with the
//! given \c key does not yet exist; otherwise an \c error_message is printed
//! and the
//! program is aborted.
template <typename Map, typename Key, typename Value>
void insertIfKeyUniqueElseError(Map& map, Key const& key, Value&& value,
                                std::string const& error_message)
{
    auto const inserted = map.emplace(key, std::forward<Value>(value));
    if (!inserted.second)
    {  // insertion failed, i.e., key already exists
        OGS_FATAL("{} Key `{}' already exists.", error_message, key);
    }
}

//! Returns the value of \c key from the given \c map if such an entry exists;
//! otherwise an \c error_message is printed and the program is aborted.
//! Cf. also the const overload below.
template <typename Map, typename Key>
OGS_NO_DANGLING typename Map::mapped_type& getOrError(
    Map& map, Key const& key, std::string const& error_message)
{
    auto it = map.find(key);
    if (it == map.end())
    {
        if constexpr (std::is_convertible<Key, std::string>::value)
        {
            OGS_FATAL("{:s} Key `{:s}' does not exist.", error_message, key);
        }
        else
        {
            OGS_FATAL("{:s} Key `{:s}' does not exist.", error_message,
                      std::to_string(key));
        }
    }

    return it->second;
}
//! \overload
template <typename Map, typename Key>
OGS_NO_DANGLING typename Map::mapped_type const& getOrError(
    Map const& map, Key const& key, std::string const& error_message)
{
    auto it = map.find(key);
    if (it == map.end())
    {
        if constexpr (std::is_convertible<Key, std::string>::value)
        {
            OGS_FATAL("{:s} Key `{:s}' does not exist.", error_message, key);
        }
        else
        {
            OGS_FATAL("{:s} Key `{:s}' does not exist.", error_message,
                      std::to_string(key));
        }
    }

    return it->second;
}

//! Returns the value of from the given \c container if such an entry fulfilling
//! the \c predicate exists;
//! otherwise an \c error_message is printed and the program is aborted.
template <typename Container, typename Predicate>
OGS_NO_DANGLING typename Container::value_type const& getIfOrError(
    Container const& container,
    Predicate&& predicate,
    std::string const& error_message)
{
    auto it = std::find_if(begin(container), end(container), predicate);
    if (it == end(container))
    {
        OGS_FATAL("Could not find element matching the predicate: {:s}",
                  error_message);
    }
    return *it;
}

/// Make the entries of the std::vector \c v unique. The remaining entries will
/// be sorted.
template <typename T>
void makeVectorUnique(std::vector<T>& v)
{
    std::sort(v.begin(), v.end());
    auto it = std::unique(v.begin(), v.end());
    v.erase(it, v.end());
}

/// Make the entries of the std::vector \c v unique using the given binary
/// function. The remaining entries will be sorted.
template <typename T, class Compare>
void makeVectorUnique(std::vector<T>& v, Compare comp)
{
    std::sort(v.begin(), v.end(), comp);
    auto it = std::unique(v.begin(), v.end());
    v.erase(it, v.end());
}

/**
 *  Reorder a vector by a given index vector.
 *
 *  Note: It is good enough in performance for medium size vectors.
 */
template <typename ValueType, typename IndexType>
void reorderVector(std::vector<ValueType>& v,
                   std::vector<IndexType> const& order)
{
    std::vector<ValueType> temp_v(v.size());
    temp_v.swap(v);

    for (std::size_t i = 0; i < order.size(); i++)
    {
        std::swap(v[i], temp_v[order[i]]);
    }
}

template <typename Container>
void uniquePushBack(Container& container,
                    typename Container::value_type const& element)
{
    if (std::find(container.begin(), container.end(), element) ==
        container.end())
    {
        container.push_back(element);
    }
}

template <typename Container>
std::optional<typename Container::value_type> findFirstNotEqualElement(
    Container const& container, typename Container::value_type const& element)
{
    auto const it =
        std::find_if_not(container.begin(), container.end(),
                         [&element](typename Container::value_type const& e)
                         { return e == element; });
    return it == container.end() ? std::nullopt : std::make_optional(*it);
}

/// Returns the index of first element in container or, if the element is not
/// found a std::size_t maximum value.
///
/// The maximum value of std::size_t is chosen, because such an index cannot
/// exist in a container; the maximum index is std::size_t::max-1.
template <typename Container>
std::size_t findIndex(Container const& container,
                      typename Container::value_type const& element)
{
    auto const it = std::find(container.begin(), container.end(), element);
    if (it == container.end())
    {
        return std::numeric_limits<std::size_t>::max();
    }
    return std::distance(container.begin(), it);
}

/** Function to destruct objects stored in a container as pointers. */
template <typename T>
void cleanupVectorElements(std::vector<T*>& items)
{
    for (auto item : items)
    {
        delete item;
    }
    items.clear();
}

/** Util function to cleanup the memory of multiple containers containing
 * pointers to objects. Sometimes, there are dependencies between the pointer
 * items in the containers. For instance, a GeoLib::Polyline or a
 * GeoLib::Surface depends on the GeoLib::Point pointers stored in a
 * std::vector<GeoLib::Point*>. Then, the dependent items have to cleaned up
 * before the GeoLib::Point objects are deleted. A similar relation exists
 * between MeshLib::Element objects and MeshLib::Node objects.*/
template <typename T1, typename... Args>
void cleanupVectorElements(std::vector<T1*>& dependent_items, Args&&... args)
{
    cleanupVectorElements(dependent_items);
    cleanupVectorElements(std::forward<Args>(args)...);
}

/// Converts range of sizes to a vector of offsets. First offset is 0 and the
/// resulting vector size is the size of the input range plus one.
template <ranges::range R>
    requires std::is_integral_v<ranges::range_value_t<R>>
std::vector<ranges::range_value_t<R>> sizesToOffsets(R const& sizes)
{
    return ranges::views::concat(
               ranges::views::single(ranges::range_value_t<R>{0}),
               ranges::views::partial_sum(sizes)) |
           ranges::to<std::vector<ranges::range_value_t<R>>>();
}

/// Checks if any of the elements in the given list is true.
template <typename List>
constexpr bool any_of(List const& values)
{
    // std::any_of is not constexpr enough in some STLs
    for (auto& value : values)
    {
        if (static_cast<bool>(value))
        {
            return true;
        }
    }

    return false;
}

/// Checks if all of the elements in the given list are true.
template <typename List>
constexpr bool all_of(List const& values)
{
    // std::all_of is not constexpr enough in some STLs
    for (auto& value : values)
    {
        if (!static_cast<bool>(value))
        {
            return false;
        }
    }

    return true;
}

/// Checks if none of the elements in the given list are true.
template <typename List>
constexpr bool none_of(List const& values)
{
    return !any_of(values);
}

/// A utility to combine multiple lambda functions into a single overloaded
/// function object.
/// Can be used with `std::visit` or similar functions requiring a callable that
/// handles multiple types.
template <class... Ts>
struct Overloaded : Ts...
{
    using Ts::operator()...;
};
#if defined(__clang__)
#if (__clang_major__ <= 16)
/// Explicit deduction guide needed for clang <= 16.
template <class... Ts>
Overloaded(Ts...) -> Overloaded<Ts...>;
#endif
#endif

}  // namespace BaseLib
