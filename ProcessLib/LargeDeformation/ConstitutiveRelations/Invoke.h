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

namespace ProcessLib::LargeDeformation
{
namespace detail
{
template <typename Result, typename Class, typename... Args>
constexpr bool areEvalArgumentTypesUnique(Result (Class::*)(Args...))
{
    using namespace boost::mp11;
    return mp_is_set<
        mp_transform<std::remove_cvref_t, mp_list<Args...>>>::value;
}

template <typename Result, typename Class, typename... Args>
constexpr bool areEvalArgumentTypesUnique(Result (Class::*)(Args...) const)
{
    using namespace boost::mp11;
    return mp_is_set<
        mp_transform<std::remove_cvref_t, mp_list<Args...>>>::value;
}
}  // namespace detail

/// Checks whether the argument types of the eval() method of the given type T
/// are unique.
///
/// Argument types differing only in constness, reference or volatility are
/// considered equal.
template <typename T>
constexpr bool areEvalArgumentTypesUnique()
{
    return detail::areEvalArgumentTypesUnique(&T::eval);
}

/// Statically asserts that the argument types of the passed Model's eval()
/// method are unique.
template <typename Model>
constexpr void assertEvalArgsUnique(Model const&)
{
    static_assert(areEvalArgumentTypesUnique<std::remove_cvref_t<Model>>());
}

}  // namespace ProcessLib::LargeDeformation
