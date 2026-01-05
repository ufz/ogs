// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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
