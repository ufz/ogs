/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <type_traits>

// ToDo (tm) remove with c++23

namespace BaseLib
{
/**
 *  \brief Converts an enumeration to its underlying type.
 *  \param e Enumeration value to convert
 *  \return The integer value of the underlying type of Enum, converted from e.
 *
 * https://en.cppreference.com/w/cpp/utility/to_underlying
 * Implementation from Scott Meyers : Modern effective c++ , Item 10 Scoped
 * enums
 */

template <typename E>
constexpr auto to_underlying(E e) noexcept
{
    return static_cast<std::underlying_type_t<E>>(e);
}
}  // namespace BaseLib