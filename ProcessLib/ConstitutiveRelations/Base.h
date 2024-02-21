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

#include "ParameterLib/SpatialPosition.h"

namespace ProcessLib::ConstitutiveRelations
{
/// Represents a previous state of type T.
template <typename T>
struct PrevState
{
    PrevState() = default;
    explicit PrevState(T const& t) : t{t} {}
    explicit PrevState(T&& t) : t{std::move(t)} {}

    PrevState<T>& operator=(T const& u)
    {
        t = u;
        return *this;
    }

    PrevState<T>& operator=(T&& u)
    {
        t = std::move(u);
        return *this;
    }

    T& operator*() { return t; }
    T const& operator*() const { return t; }

    T* operator->() { return &t; }
    T const* operator->() const { return &t; }

private:
    T t;
};

//! Applies PrevState to a tuple of constitutive data.
template <typename Tuple>
using PrevStateOf = boost::mp11::mp_transform<PrevState, Tuple>;

namespace detail
{
template <typename... Ts, std::size_t... Idcs>
void assign(std::tuple<PrevState<Ts>...>& prev_states,
            std::tuple<Ts...> const& current_states,
            std::index_sequence<Idcs...>)
{
    ((std::get<Idcs>(prev_states) = std::get<Idcs>(current_states)), ...);
}
}  // namespace detail

//! Assigns a tuple of current states to a tuple of previous states.
template <typename... Ts>
void assign(std::tuple<PrevState<Ts>...>& prev_states,
            std::tuple<Ts...> const& current_states)
{
    detail::assign(prev_states, current_states,
                   std::make_index_sequence<sizeof...(Ts)>{});
}

struct SpaceTimeData
{
    ParameterLib::SpatialPosition x;
    double t;
    double dt;
};

/// Convenience alias for not a number.
static constexpr double nan = std::numeric_limits<double>::quiet_NaN();

}  // namespace ProcessLib::ConstitutiveRelations
