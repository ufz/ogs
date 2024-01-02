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

#include <string>
#include <tuple>

namespace ProcessLib::Reflection
{
/**
 * Provides access to a single member of \c Class via an accessor function and
 * provides a name that can be used, e.g. to identify that member during I/O
 * operations.
 *
 * \c Accessor is a function taking an instance of \c Class and returning a
 * reference to the member.
 */
template <typename Class, typename Accessor>
struct ReflectionData
{
    static_assert(std::is_same_v<Accessor, std::remove_cvref_t<Accessor>>);
    static_assert(std::is_invocable_v<Accessor, Class&>);
    static_assert(std::is_invocable_v<Accessor, Class const&>);
    static_assert(
        std::is_lvalue_reference_v<std::invoke_result_t<Accessor, Class&>>);
    static_assert(std::is_lvalue_reference_v<
                  std::invoke_result_t<Accessor, Class const&>>);

    ReflectionData(std::string name, Accessor&& accessor)
        : name(std::move(name)), accessor(std::move(accessor))
    {
    }

    explicit ReflectionData(Accessor&& accessor) : accessor(std::move(accessor))
    {
    }

    std::string name;
    Accessor accessor;
};

template <typename Class, typename Accessor>
auto makeReflectionData(Accessor&& accessor)
{
    return ReflectionData<Class, std::remove_cvref_t<Accessor>>{
        std::forward<Accessor>(accessor)};
}

template <typename Class, typename Accessor>
auto makeReflectionData(std::string name, Accessor&& accessor)
{
    return ReflectionData<Class, std::remove_cvref_t<Accessor>>{
        std::move(name), std::forward<Accessor>(accessor)};
}

template <typename Class, typename Member>
auto makeReflectionData(Member Class::*member)
{
    return makeReflectionData<Class>([member](auto& obj) -> auto&
                                     { return obj.*member; });
}

template <typename Class, typename Member>
auto makeReflectionData(std::string name, Member Class::*member)
{
    return makeReflectionData<Class>(
        std::move(name), [member](auto& obj) -> auto& { return obj.*member; });
}

template <typename Class, typename Member>
auto reflectWithName(std::string name, Member Class::*member)
{
    return std::tuple{makeReflectionData(std::move(name), member)};
}

template <typename Class, typename... Accessors>
auto reflectWithoutName(Accessors&&... accessors)
{
    return std::tuple{
        makeReflectionData<Class>(std::forward<Accessors>(accessors))...};
}

template <typename Class, typename... Members>
auto reflectWithoutName(Members Class::*... members)
{
    return std::tuple{makeReflectionData<Class>(members)...};
}
}  // namespace ProcessLib::Reflection
