/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <tuple>
#include <utility>

#include "Apply.h"

namespace ProcessLib::Graph
{
namespace detail
{
template <typename T>
concept HasCreate = requires
{
    T::create;
};

template <typename Model, typename TupleOfArgs>
Model constructModel(TupleOfArgs& args)
{
    if constexpr (HasCreate<Model>)
    {
        return applyImpl(&Model::create, args);
    }
    else if constexpr (std::is_default_constructible_v<Model>)
    {
        return Model{};
    }
    else
    {
        static_assert(std::is_default_constructible_v<
                          Model> /* This is definitely false, here. */,
                      "The model to be constructed has neither a static "
                      "create() function nor is it default constructible.");
    }
}
template <template <typename...> typename Tuple,
          typename... Models,
          typename TupleOfArgs>
Tuple<Models...> constructModels(std::type_identity<Tuple<Models...>>,
                                 TupleOfArgs&& args)
{
    return Tuple{constructModel<Models>(args)...};
}
}  // namespace detail

/**
 * Constructs a tuple of models.
 *
 * Each model in the tuple is either
 *
 * 1. constructed by calling a static create() method of the model's class,
 * 2. or by default construction if the former does not exist.
 *
 * In case (1.) the arguments are passed via their data types from the passed
 * \c args to the create() method.
 */
template <typename TupleOfModels, typename... Args>
TupleOfModels constructModels(Args&&... args)
{
    return detail::constructModels(
        std::type_identity<TupleOfModels>{},
        std::forward_as_tuple(std::forward<Args>(args)...));
}
}  // namespace ProcessLib::Graph
