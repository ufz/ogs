/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <boost/mp11.hpp>

namespace ProcessLib::Graph
{
namespace detail
{
template <typename T, typename Tuple, typename... Tuples>
auto& getImpl(Tuple& t, Tuples&... ts)
{
    using namespace boost::mp11;

    using PlainTuple = std::remove_cv_t<Tuple>;
    using TupleOfPlainTypes = mp_transform<std::remove_cvref_t, PlainTuple>;
    using Index = mp_find<TupleOfPlainTypes, T>;

    if constexpr (Index::value != mp_size<PlainTuple>::value)
    {
        return std::get<Index::value>(t);
    }
    else
    {
        return getImpl<T>(ts...);
    }
}

template <typename... Tuples>
struct GetFlattenedTupleTypes
{
    // This check is necessary only to ensure that mp_flatten works as expected.
    static_assert(boost::mp11::mp_similar<std::tuple<>, Tuples...>::value,
                  "All template arguments must be std::tuple's.");

    using type = boost::mp11::mp_flatten<std::tuple<Tuples...>>;
};

// Alias to be used together with static_assert to improve possible compile
// error messages.
template <typename List, typename Elem>
constexpr bool mp_contains_v = boost::mp11::mp_contains<List, Elem>::value;

// Alias to be used together with static_assert to improve possible compile
// error messages
template <typename Set>
constexpr bool mp_is_set_v = boost::mp11::mp_is_set<Set>::value;
}  // namespace detail

/// Type-based access of an element of any of the passed tuples.
///
/// This function does essentially the same as
/// <code>std::get<T>(some_tuple)</code>, but for any number of passed tuples.
///
/// The type \c T must be present in the \c Tuples's member types exactly  once.
/// The passed \c Tuples's member types might be cvref qualified, but \c T must
/// not.
template <typename T, typename... Tuples>
auto& get(Tuples&... ts)
{
    using namespace boost::mp11;

    static_assert(std::is_same_v<T, std::remove_cvref_t<T>>,
                  "The passed type T must not be cvref qualified.");

    using FlattenedTuple = typename detail::GetFlattenedTupleTypes<
        std::remove_cv_t<Tuples>...>::type;
    using FlattenedTupleOfPlainTypes =
        mp_transform<std::remove_cvref_t, FlattenedTuple>;

    static_assert(
        detail::mp_is_set_v<FlattenedTupleOfPlainTypes>,
        "The types of all elements of all passed tuples must be unique.");

    static_assert(detail::mp_contains_v<FlattenedTupleOfPlainTypes, T>,
                  "Type T must be inside any of the passed tuples.");

    return detail::getImpl<T>(ts...);
}

}  // namespace ProcessLib::Graph
