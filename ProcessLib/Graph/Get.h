/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include "boost/mp11.hpp"

namespace ProcessLib::Graph
{
namespace detail
{
template <typename T, typename Tuple, typename... Tuples>
auto& getImpl(Tuple& t, Tuples&... ts)
{
    using namespace boost::mp11;

    if constexpr (mp_contains<std::remove_cv_t<Tuple>, T>::value)
    {
        return std::get<T>(t);
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
}  // namespace detail

/// Type-based access of an element of any of the passed tuples.
///
/// This function does essentially the same as
/// <code>std::get<T>(some_tuple)</code>, but for any number of passed tuples.
///
/// The type \code T must be present in the \code Tuples's types exactly once.
template <typename T, typename... Tuples>
auto& get(Tuples&... ts)
{
    using namespace boost::mp11;

    using FlattenedTuple = typename detail::GetFlattenedTupleTypes<
        std::remove_cvref_t<Tuples>...>::type;

    static_assert(
        mp_is_set<FlattenedTuple>::value,
        "The types of all elements of all passed tuples must be unique.");

    static_assert(mp_contains<FlattenedTuple, T>::value,
                  "Type T must be inside any of the passed tuples.");

    return detail::getImpl<T>(ts...);
}

}  // namespace ProcessLib::Graph
