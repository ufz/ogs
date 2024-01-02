/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <concepts>
#include <utility>

namespace BaseLib
{
//! A "strong" version of the underlying type T.
//!
//! Strong types with same T but different tag are not implicitly convertible to
//! each other, and are not implicitly convertible to/from T. Based on these
//! properties strong types can serve as an ingredient for safer APIs.
template <typename T, typename Tag>
struct StrongType
{
    static_assert(
        !std::is_reference_v<T>,
        "The underlying type must not be a reference. This struct is not "
        "designed for references. If its design should be extended to also "
        "cover references make sure to unit test it throroughly. In this case "
        "also consider the implications for the overall usage of StrongType in "
        "the code, e.g., StrongType<int, struct Tag> and StrongType<int&, "
        "struct Tag> are not related to each other by any inheritance "
        "relation, but are completely different types, and in particular "
        "cannot be passed to functions expecting a value/reference of the "
        "other type. Therefore, allowing references as underlying types might "
        "have adverse effect on \"large scale\" API design.");

    static_assert(
        std::default_initializable<T>,
        "The underlying type must be default initializable. This struct is "
        "designed for \"general purpose\" underlying types such as numbers or "
        "matrices, which are usually default constructible. If you need "
        "support for non-default-constuctible underlying types, make sure you "
        "test your extension of the current implementation thoroughly.");

    constexpr StrongType() noexcept(std::is_nothrow_default_constructible_v<T>)
        : value_{}  // "zero initialization"
    {
    }

    constexpr explicit StrongType(T const& value) noexcept(
        std::is_nothrow_copy_constructible_v<T>)
        : value_{value}
    {
    }

    constexpr explicit StrongType(T&& value) noexcept(
        std::is_nothrow_move_constructible_v<T>)
        : value_{std::move(value)}
    {
    }

    //! Value access. Alternative to operator*, operator() might look clearer in
    //! mathematical expressions.
    [[nodiscard]] constexpr T& operator()() noexcept { return this->value_; }
    [[nodiscard]] constexpr T const& operator()() const noexcept
    {
        return value_;
    }

    //! Value access. Alternative to operator(). Introduced for symmetry reasons
    //! as a companion to operator->.
    [[nodiscard]] constexpr T& operator*() noexcept { return this->value_; }
    [[nodiscard]] constexpr T const& operator*() const noexcept
    {
        return value_;
    }

    //! Value access.
    [[nodiscard]] constexpr T* operator->() noexcept { return &this->value_; }
    [[nodiscard]] constexpr T const* operator->() const noexcept
    {
        return &value_;
    }

private:
    T value_;
};
}  // namespace BaseLib
