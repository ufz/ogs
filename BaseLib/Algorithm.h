/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <algorithm>
#include <cassert>
#include <typeindex>
#include <typeinfo>
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
std::vector<T> excludeObjectCopy(std::vector<T> const& src_vec,
    std::vector<std::size_t> const& exclude_positions)
{
    std::vector<T> dest_vec;
    if (exclude_positions.empty()) {
        dest_vec = src_vec;
        return dest_vec;
    }

    assert (exclude_positions.back() < src_vec.size());

    std::copy_n(src_vec.cbegin(), exclude_positions[0], std::back_inserter(dest_vec));
    for (std::size_t i=1; i<exclude_positions.size(); ++i) {
        std::copy_n(
            src_vec.cbegin()+exclude_positions[i-1]+1,
            exclude_positions[i]-(exclude_positions[i-1]+1),
            std::back_inserter(dest_vec)
        );
    }
    std::copy(src_vec.cbegin()+exclude_positions.back()+1,
        src_vec.cend(), std::back_inserter(dest_vec));

    return dest_vec;
}

template <typename T>
void excludeObjectCopy(std::vector<T> const& src_vec,
    std::vector<std::size_t> const& exclude_positions,
    std::vector<T> & dest_vec)
{
    dest_vec = excludeObjectCopy(src_vec, exclude_positions);
}

template <typename InputIt, typename Predicate>
typename std::iterator_traits<InputIt>::reference findElementOrError(
    InputIt begin, InputIt end, Predicate predicate,
    std::string const& error = "")
{
    auto it = std::find_if(begin, end, predicate);
    if (it == end)
    {
        OGS_FATAL("Element not found in the input range; %s", error.c_str());
    }
    return *it;
}

//! Inserts the given \c key with the given \c value into the \c map if an entry
//! with the
//! given \c key does not yet exist; otherwise an \c error_message is printed
//! and the
//! program is aborted.
//! Note: The type of \c key must be std::type_index.
template <typename Map, typename Key, typename Value>
void insertIfTypeIndexKeyUniqueElseError(Map& map, Key const& key,
                                         Value&& value,
                                         std::string const& error_message)
{
    auto const inserted = map.emplace(key, std::forward<Value>(value));
    if (!inserted.second)
    {  // insertion failed, i.e., key already exists
        OGS_FATAL("%s Key `%s' already exists.", error_message.c_str(),
                  tostring(key.hash_code()).c_str());
    }
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
        OGS_FATAL("%s Key `%s' already exists.", error_message.c_str(),
                  tostring(key).c_str());
    }
}

//! Inserts the given \c key with the given \c value into the \c map if neither
//! an entry
//! with the given \c key nor an entry with the given \c value already exists;
//! otherwise an \c error_message is printed and the program is aborted.
template <typename Map, typename Key, typename Value>
void insertIfKeyValueUniqueElseError(Map& map, Key const& key, Value&& value,
                                     std::string const& error_message)
{
    auto value_compare = [&value](typename Map::value_type const& elem) {
        return value == elem.second;
    };

    if (std::find_if(map.cbegin(), map.cend(), value_compare) != map.cend())
    {
        OGS_FATAL("%s Value `%s' already exists.", error_message.c_str(),
                  tostring(value).c_str());
    }

    auto const inserted = map.emplace(key, std::forward<Value>(value));
    if (!inserted.second)
    {  // insertion failed, i.e., key already exists
        OGS_FATAL("%s Key `%s' already exists.", error_message.c_str(),
                  tostring(key).c_str());
    }
}

//! Returns the value of \c key from the given \c map if such an entry exists;
//! otherwise an \c error_message is printed and the program is aborted.
//! Cf. also the const overload below.
//! \remark Use as: \code{.cpp} get_or_error<Value>(some_map, some_key, "error
//! message") \endcode
template <typename Map, typename Key>
typename Map::mapped_type& getOrError(Map& map, Key const& key,
                                      std::string const& error_message)
{
    auto it = map.find(key);
    if (it == map.end())
    {
        OGS_FATAL("%s Key `%s' does not exist.", error_message.c_str(),
                  tostring(key).c_str());
    }

    return it->second;
}
//! \overload
template <typename Map, typename Key>
typename Map::mapped_type const& getOrError(Map const& map, Key const& key,
                                            std::string const& error_message)
{
    auto it = map.find(key);
    if (it == map.end())
    {
        OGS_FATAL("%s Key `%s' does not exist.", error_message.c_str(),
                  tostring(key).c_str());
    }

    return it->second;
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

    for (std::size_t i=0; i<order.size(); i++)
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
        container.push_back(element);
}
}  // namespace BaseLib
