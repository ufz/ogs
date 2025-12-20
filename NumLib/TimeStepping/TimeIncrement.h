// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <algorithm>
#include <cmath>
#include <compare>
#include <iomanip>
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
