/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <compare>
#include <limits>

#include "MathLib/KahanSum.h"

namespace NumLib
{
struct TimeIncrement;
}

template <>
struct fmt::formatter<NumLib::TimeIncrement> : fmt::ostream_formatter
{
};

namespace NumLib
{
struct TimeIncrement
{
    constexpr explicit TimeIncrement(double const dt) : value_{dt} {}

    constexpr double operator()() const { return value_; }

    friend inline std::ostream& operator<<(std::ostream& os,
                                           TimeIncrement const& dt)
    {
        auto const precision = os.precision();
        return os << std::setprecision(
                         std::numeric_limits<double>::max_digits10)
                  << dt.value_ << std::setprecision(precision);
    }

private:
    double value_;
};

}  // namespace NumLib
