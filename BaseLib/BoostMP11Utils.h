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

namespace boost::mp11
{
// Alias to be used together with static_assert to improve possible compile
// error messages.
template <typename List, typename Elem>
constexpr bool mp_contains_v = mp_contains<List, Elem>::value;

// Alias to be used together with static_assert to improve possible compile
// error messages
template <typename Set>
constexpr bool mp_is_list_v = mp_is_list<Set>::value;

// Alias to be used together with static_assert to improve possible compile
// error messages
template <typename Set>
constexpr bool mp_is_set_v = mp_is_set<Set>::value;
}  // namespace boost::mp11
