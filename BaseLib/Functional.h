/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <functional>
#include "baselib_export.h"

namespace BaseLib
{
namespace detail
{
//! Helper struct used to make std::placeholders::_1, ... accessible via a
//! compile-time computed index (the template parameter).
template <int>
struct IndexedPlacedPlaceholder;

//! Creates specializations of IndexedPlacedPlaceholder.
//! \param INDEX the integer value which is specialized
//! \param INDEX_P_1 "index plus one"; if INDEX_P_1 equals 1, then the member
//! value will be std::placeholders::_1, etc.
#define SPECIALIZE_INDEXEDPLACEHOLDER(INDEX, INDEX_P_1)               \
    template <>                                                       \
    struct IndexedPlacedPlaceholder<(INDEX)> {                        \
        static BASELIB_EXPORT const decltype(std::placeholders::_##INDEX_P_1) value; \
    }

// Create specializations up to the tenth placeholder
SPECIALIZE_INDEXEDPLACEHOLDER(0, 1);
SPECIALIZE_INDEXEDPLACEHOLDER(1, 2);
SPECIALIZE_INDEXEDPLACEHOLDER(2, 3);
SPECIALIZE_INDEXEDPLACEHOLDER(3, 4);
SPECIALIZE_INDEXEDPLACEHOLDER(4, 5);
SPECIALIZE_INDEXEDPLACEHOLDER(5, 6);
SPECIALIZE_INDEXEDPLACEHOLDER(6, 7);
SPECIALIZE_INDEXEDPLACEHOLDER(7, 8);
SPECIALIZE_INDEXEDPLACEHOLDER(8, 9);
SPECIALIZE_INDEXEDPLACEHOLDER(9, 10);

#undef SPECIALIZE_INDEXEDPLACEHOLDER

// Note: The call sequence is easyBind() -> easyBind_inner() ->
// easyBind_innermost().

template <int... Indices, typename Object, typename ReturnType,
          typename... Args>
std::function<ReturnType(Args...)> easyBind_innermost(
    ReturnType (Object::*method)(Args...), Object& obj)
{
    // std::ref makes sure that obj is not copied.
    return std::bind(method, std::ref(obj),
                     IndexedPlacedPlaceholder<Indices>::value...);
}

template <int... Indices, typename Object, typename ReturnType,
          typename... Args>
std::function<ReturnType(Args...)> easyBind_innermost(
    ReturnType (Object::*method)(Args...) const, Object const& obj)
{
    // std::cref makes sure that obj is not copied.
    return std::bind(method, std::cref(obj),
                     IndexedPlacedPlaceholder<Indices>::value...);
}

template <int... Indices, typename Object, typename ReturnType,
          typename... Args>
std::function<ReturnType(Args...)> easyBind_innermost(
    ReturnType (Object::*method)(Args...) const, Object& obj)
{
    // std::cref makes sure that obj is not copied.
    return std::bind(method, std::cref(obj),
                     IndexedPlacedPlaceholder<Indices>::value...);
}

template <int... Indices, typename Object, typename MethodClass,
          typename ReturnType, typename... Args>
std::function<ReturnType(Args...)> easyBind_innermost(
    ReturnType (MethodClass::*method)(Args...), Object&& obj)
{
    return std::bind(method, std::forward<Object>(obj),
                     IndexedPlacedPlaceholder<Indices>::value...);
}

template <int... Indices, typename Object, typename MethodClass,
          typename ReturnType, typename... Args>
std::function<ReturnType(Args...)> easyBind_innermost(
    ReturnType (MethodClass::*method)(Args...) const, Object&& obj)
{
    return std::bind(method, std::forward<Object>(obj),
                     IndexedPlacedPlaceholder<Indices>::value...);
}

template <int... Indices, typename Object, typename MethodClass,
          typename ReturnType, typename... Args>
std::function<ReturnType(Args...)> easyBind_inner(
    ReturnType (MethodClass::*method)(Args...), Object&& obj,
    std::integer_sequence<int, Indices...>)
{
    return easyBind_innermost<Indices...>(method, std::forward<Object>(obj));
}

template <int... Indices, typename Object, typename MethodClass,
          typename ReturnType, typename... Args>
std::function<ReturnType(Args...)> easyBind_inner(
    ReturnType (MethodClass::*method)(Args...) const, Object&& obj,
    std::integer_sequence<int, Indices...>)
{
    return easyBind_innermost<Indices...>(method, std::forward<Object>(obj));
}

/*! Deduces the signature of the call operator of class \c T.
 *
 * The matching type of std::function is provided as the member type
 * \c FunctionType.
 *
 * \see http://stackoverflow.com/a/7943765
 */
template <typename T>
struct FunctionTraits
    : public FunctionTraits<decltype(&std::decay<T>::type::operator())> {
};

template <typename Object, typename ReturnType, typename... Args>
struct FunctionTraits<ReturnType (Object::*)(Args...)> {
    using FunctionType = std::function<ReturnType(Args...)>;
};

template <typename Object, typename ReturnType, typename... Args>
struct FunctionTraits<ReturnType (Object::*)(Args...) const> {
    using FunctionType = std::function<ReturnType(Args...)>;
};

}  // namespace detail

/*! Convenience wrapper for std::bind().
 *
 * This function binds the member function pointer \c member of class \c Object
 * to the instance \c obj of this class and wraps the result in a std::function
 * with matching signature.
 *
 * The result of this function can be used, e.g., to deduce the signature of the
 * \c method (which is not possible with the result of std::bind).
 *
 * Example:
 * \code{.cpp}
 * using std::placeholders;
 * Object some_object;
 *
 * auto f_bind = std::function<ReturnType(Arg1, Arg2, Arg3>(
 *                  std::bind(&Object::methodWithThreeArguments,
 *                      std::ref(some_object), _1, _2, _3));
 *
 * auto f_easy = easyBind(&Object::methodWithThreeArguments, some_object);
 * \endcode
 *
 * In the example the expressions creating \c f_bind and \c f_easy are
 * equivalent.
 *
 * \note
 * There is one difference between the behaviour of std::bind and the one of
 * easyBind: In easyBind \c obj is never copied, instead it will be referenced.
 * This is in contrast to the behaviour of std::bind, and has been chosen in
 * order to prevent accidental copies.
 */
template <typename Object, typename MethodClass, typename ReturnType,
          typename... Args>
typename std::enable_if<
    std::is_same<MethodClass,
                 /* Note: All of remove_cv, remove_pointer and decay is
                  * necessary, e.g. if method is a member function pointer. */
                 typename std::remove_cv<typename std::remove_pointer<
                     typename std::decay<Object>::type>::type>::type>::value,
    std::function<ReturnType(Args...)>>::type
easyBind(ReturnType (MethodClass::*method)(Args...), Object&& obj)
{
    return detail::easyBind_inner(
        method, std::forward<Object>(obj),
        std::make_integer_sequence<int, sizeof...(Args)>{});
}

//! \overload
template <typename Object, typename MethodClass, typename ReturnType,
          typename... Args>
typename std::enable_if<
    std::is_same<MethodClass,
                 typename std::remove_cv<typename std::remove_pointer<
                     typename std::decay<Object>::type>::type>::type>::value,
    std::function<ReturnType(Args...)>>::type
easyBind(ReturnType (MethodClass::*method)(Args...) const, Object&& obj)
{
    return detail::easyBind_inner(
        method, std::forward<Object>(obj),
        std::make_integer_sequence<int, sizeof...(Args)>{});
}

//! Wraps a callable object in a std::function.
//!
//! This method is provided for convenience since it automatically deduces the
//! correct type of std::function.
template <typename Object>
typename detail::FunctionTraits<Object>::FunctionType easyBind(Object&& obj)
{
    return BaseLib::easyBind(&std::decay<Object>::type::operator(),
                             std::forward<Object>(obj));
}

}  // namespace BaseLib
