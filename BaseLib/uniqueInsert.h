/**
 * \file
 * \author Thomas Fischer
 * \date   2011-02-23
 * \brief  Definition of the uniquePushBack function.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <algorithm>
#include <typeinfo>
#include <typeindex>
#include "Error.h"

namespace BaseLib
{
template <typename Container>
void uniquePushBack(Container& container,
                    typename Container::value_type const& element)
{
    if (std::find(container.begin(), container.end(), element) ==
        container.end())
        container.push_back(element);
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

}  // end namespace BaseLib
