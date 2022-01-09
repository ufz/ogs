/**
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <algorithm>
#include <array>
#include <utility>

#include "Algorithm.h"

namespace BaseLib::TMP
{
// Concat ----------------------------------------------------------------------
/**
 * Concatenates two lists of types.
 */
template <typename List1, typename List2>
struct Concat;

template <template <typename... /*Types*/> typename List, typename... Types1,
          typename... Types2>
struct Concat<List<Types1...>, List<Types2...>>
{
    using type = List<Types1..., Types2...>;
};

/// \copydoc BaseLib::TMP::Concat
template <typename List1, typename List2>
using Concat_t = typename Concat<List1, List2>::type;

// filter- ---------------------------------------------------------------------
/// A list of types.
template <typename...>
struct List
{
};

namespace detail
{
template <typename Pred, template <typename...> typename SomeListOfTypes>
constexpr List<> filter(Pred, SomeListOfTypes<>*)
{
    return {};
}

template <typename Pred, template <typename...> typename SomeListOfTypes,
          typename Head, typename... Tail>
constexpr decltype(auto) filter(Pred pred, SomeListOfTypes<Head, Tail...>*)
{
    constexpr auto tail_arg = static_cast<List<Tail...>*>(nullptr);
    using TailFiltered = decltype(filter(pred, tail_arg));

    using Result =
        std::conditional_t<pred((Head*)nullptr),
                           Concat_t<List<Head>, TailFiltered>, TailFiltered>;

    return Result{};
}
}  // namespace detail

/**
 * Filters the given list of types via the given predicate.
 *
 * Keeps all elements for which \c pred is true.
 *
 * The function signature of \c pred  should be <tt>bool pred(Type*)</tt>, where
 * \c Type is a member of the given list of types.
 */
template <typename List, typename Pred>
decltype(auto) filter(Pred pred)
{
    return detail::filter(pred, static_cast<List*>(nullptr));
}

// Map -------------------------------------------------------------------------
/**
 * Turns a list of types into a list of different types via the provided type
 * mapping.
 */
template <template <typename /*FromType*/> typename MapFromTypeToType,
          typename List>
struct Map;

template <template <typename /*FromType*/> typename MapFromTypeToType,
          template <typename... /*Types*/> typename List, typename... Types>
struct Map<MapFromTypeToType, List<Types...>>
{
    using type = List<MapFromTypeToType<Types>...>;
};

/// \copydoc BaseLib::TMP::Map
template <template <typename /*FromType*/> typename MapFromTypeToType,
          typename List>
using Map_t = typename Map<MapFromTypeToType, List>::type;

// map_to_array ----------------------------------------------------------------
namespace detail
{
template <template <typename... /*Types*/> typename List, typename Function,
          typename Head, typename... Tail>
constexpr decltype(auto) map_to_array(Function&& f, List<Head, Tail...>*)
{
    using CommonType =
        std::common_type_t<std::invoke_result_t<Function, Head*>,
                           std::invoke_result_t<Function, Tail*>...>;

    return std::array<CommonType, 1 + sizeof...(Tail)>{f((Head*)nullptr),
                                                       f((Tail*)nullptr)...};
}

template <template <typename... /*Types*/> typename List, typename Function>
constexpr std::array<std::nullptr_t, 0> map_to_array(Function&& /*f*/, List<>*)
{
    return {};
}
}  // namespace detail

/**
 * Applies the given function to a list of types returning the array of results.
 *
 * The function signature should be <tt>ReturnType f(Type*)</tt>, where \c Type
 * is a member of the given list of types and \c ReturnType is the same for all
 * overloads of \c f.
 */
template <typename List, typename Function>
constexpr decltype(auto) map_to_array(Function&& f)
{
    return detail::map_to_array(std::forward<Function>(f), (List*)nullptr);
}

// foreach ---------------------------------------------------------------------
namespace detail
{
template <template <typename... /*Types*/> typename List, typename Function,
          typename... Types>
void foreach (Function&& f, List<Types...>*)
{
    (..., f((Types*)nullptr));
}
}  // namespace detail

/**
 * Applies the given function to a list of types.
 *
 * The function signature should be <tt>void f(Type*)</tt>, where \c Type is a
 * member of the given list of types.
 */
template <typename List, typename Function>
void foreach (Function&& f)
{
    detail::foreach (std::forward<Function>(f), (List*)nullptr);
}

// contains --------------------------------------------------------------------
/// Returns if \c Type is contained in the given \c List of types.
template <typename List, typename Type>
constexpr bool contains()
{
    auto pred_is_same = []<typename OtherType>(OtherType*)
    { return std::is_same_v<Type, OtherType>; };

    return any_of(map_to_array<List>(pred_is_same));
}

}  // namespace BaseLib::TMP
