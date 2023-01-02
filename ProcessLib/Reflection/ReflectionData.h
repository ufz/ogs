/**
 * \file
 * \copyright
 * Copyright (c) 2012-2023, OpenGeoSys Community (http://www.opengeosys.org)
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
template <typename Class, typename Member>
struct ReflectionData
{
    ReflectionData(std::string name, Member Class::*field)
        : name(std::move(name)), field(field)
    {
    }

    explicit ReflectionData(Member Class::*field) : field(field) {}

    std::string name;
    Member Class::*field;
};

template <typename Class, typename Member>
std::tuple<ReflectionData<Class, Member>> reflectWithName(
    std::string const& name, Member Class::*field)
{
    return {{name, field}};
}

template <typename Class, typename... Members>
auto reflectWithoutName(Members Class::*... members)
{
    return std::tuple{ReflectionData<Class, Members>{"", members}...};
}
}  // namespace ProcessLib::Reflection
