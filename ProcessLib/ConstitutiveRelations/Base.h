// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <boost/mp11.hpp>

#include "BaseLib/StrongType.h"
#include "MathLib/KelvinVector.h"
#include "ParameterLib/SpatialPosition.h"

namespace ProcessLib::ConstitutiveRelations
{
namespace KV = MathLib::KelvinVector;

template <int DisplacementDim>
using KelvinVector = KV::KelvinVectorType<DisplacementDim>;

template <int DisplacementDim>
using KelvinMatrix = KV::KelvinMatrixType<DisplacementDim>;

template <int DisplacementDim>
using GlobalDimVector = Eigen::Vector<double, DisplacementDim>;

template <int DisplacementDim>
using GlobalDimMatrix =
    Eigen::Matrix<double, DisplacementDim, DisplacementDim, Eigen::RowMajor>;

using Temperature = BaseLib::StrongType<double, struct TemperatureTag>;

template <int DisplacementDim>
using SpecificBodyForce = BaseLib::StrongType<GlobalDimVector<DisplacementDim>,
                                              struct SpecificBodyForceTag>;

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
