// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "MaterialLib/MPL/Properties/Sigmoid.h"

#include <cmath>
#include <limits>
#include <numbers>
#include <optional>
#include <utility>

#include "BaseLib/Error.h"

namespace MaterialPropertyLib
{
// Largest exponent x for which std::exp(x) does not overflow a double.
static constexpr double max_exp_argument =
    std::numeric_limits<double>::max_exponent * std::numbers::ln2;

// value() is time-step independent; the derivatives pass this placeholder when
// they call it to obtain the current sigmoid value.
static constexpr double unused_dt = std::numeric_limits<double>::quiet_NaN();

// Shared derivative preamble: returns (f_normalized, range) at the current
// state, or nullopt when the sigmoid is saturated (f_normalized at 0 or range),
// where all derivatives vanish.
static std::optional<std::pair<double, double>> normalizedState(
    double const f, double const lower_bound, double const upper_bound)
{
    double const range = upper_bound - lower_bound;
    double const f_normalized = f - lower_bound;

    if (f_normalized <= 0.0 || f_normalized >= range)
    {
        return std::nullopt;
    }
    return std::pair{f_normalized, range};
}

Sigmoid::Sigmoid(std::string name,
                 double const steepness,
                 double const characteristic_value,
                 double const lower_bound,
                 double const upper_bound,
                 Variable const independent_variable)
    : steepness_(steepness),
      characteristic_value_(characteristic_value),
      lower_bound_(lower_bound),
      upper_bound_(upper_bound),
      independent_variable_(independent_variable)
{
    if (lower_bound_ >= upper_bound_)
    {
        OGS_FATAL(
            "In the sigmoid property '{:s}', the lower_bound value {} must be "
            "smaller than the upper_bound value {}.",
            name, lower_bound_, upper_bound_);
    }

    name_ = std::move(name);
}

PropertyDataType Sigmoid::value(VariableArray const& variable_array,
                                ParameterLib::SpatialPosition const& /*pos*/,
                                double const /*t*/,
                                double const /*dt*/) const
{
    double const X = std::get<double>(variable_array[independent_variable_]);

    // Compute sigmoid: (X_u - X_l) / (1 + exp(-k * (X - X_c))) + X_l
    double const exponent = -steepness_ * (X - characteristic_value_);

    // Exponent large positive -> sigmoid -> lower_bound; avoid exp overflow.
    if (exponent > max_exp_argument)
    {
        return lower_bound_;
    }

    double const sigmoid_value =
        (upper_bound_ - lower_bound_) / (1.0 + std::exp(exponent)) +
        lower_bound_;

    return sigmoid_value;
}

bool Sigmoid::derivativeVanishes(Variable const variable) const
{
    // Non-zero only for the derivative w.r.t. the independent variable.
    return variable != independent_variable_;
}

PropertyDataType Sigmoid::dValue(VariableArray const& variable_array,
                                 Variable const variable,
                                 ParameterLib::SpatialPosition const& pos,
                                 double const t,
                                 double const /*dt*/) const
{
    if (derivativeVanishes(variable))
    {
        return 0.0;
    }

    // Compute derivative using a numerically stable formula.
    // From dz/dX = k * z * (1-z) where z = 1 / (1 + exp(-k*(X-Xc))),
    // and f_norm = (upper - lower) * z, we get:
    // f'(X) = k * f_norm / (upper - lower) * ((upper - lower) - f_norm)
    // This avoids computing huge (1 + exp(-k*(X-Xc)))^2
    double const f = std::get<double>(value(variable_array, pos, t, unused_dt));
    auto const state = normalizedState(f, lower_bound_, upper_bound_);
    if (!state)
    {
        return 0.0;
    }
    auto const [f_normalized, range] = *state;

    return steepness_ * f_normalized * (range - f_normalized) / range;
}

PropertyDataType Sigmoid::d2Value(VariableArray const& variable_array,
                                  Variable const variable1,
                                  Variable const variable2,
                                  ParameterLib::SpatialPosition const& pos,
                                  double const t,
                                  double const /*dt*/) const
{
    // Only compute second derivative if both variables are the independent
    // variable.
    if (derivativeVanishes(variable1) || derivativeVanishes(variable2))
    {
        return 0.0;
    }

    // Compute second derivative using a numerically stable formula.
    // From f'(X) = k * f_norm * (range - f_norm) / range, we get:
    // f''(X) = k^2 * f_norm * (range - f_norm) * (range - 2*f_norm) / range^2
    // This avoids overflow/underflow issues.
    double const f = std::get<double>(value(variable_array, pos, t, unused_dt));
    auto const state = normalizedState(f, lower_bound_, upper_bound_);
    if (!state)
    {
        return 0.0;
    }
    auto const [f_normalized, range] = *state;

    return steepness_ * steepness_ * f_normalized * (range - f_normalized) *
           (range - 2.0 * f_normalized) / (range * range);
}

}  // namespace MaterialPropertyLib
