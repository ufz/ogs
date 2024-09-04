/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
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
struct Time;
}

template <>
struct fmt::formatter<NumLib::Time> : fmt::ostream_formatter
{
};

namespace NumLib
{
struct Time
{
    constexpr explicit Time(double const time) : value_{time} {}

    constexpr double operator()() const { return value_(); }

    constexpr Time& operator+=(double const increment)
    {
        value_ += increment;
        return *this;
    }

    constexpr Time operator+(double const increment) const
    {
        Time result = *this;
        result += increment;
        return result;
    }

    constexpr Time& operator-=(double const decrement)
    {
        value_ -= decrement;
        return *this;
    }

    constexpr Time operator-(double const decrement) const
    {
        Time result = *this;
        result -= decrement;
        return result;
    }

    constexpr inline std::weak_ordering operator<=>(Time const& b) const
    {
        Time const& a = *this;
        double const diff = b() - a();

        double const eps = 10 * std::numeric_limits<double>::epsilon() *
                           std::max(1., (std::abs(a()) + std::abs(b())) / 2);

        // Keep for debugging
        // DBUG("Compare {} to {}. x is {} and y is {}. Eps {}, diff {}", a, b,
        //      x, y, eps, diff);

        if (diff < -eps)
        {
            return std::weak_ordering::greater;
        }
        if (diff > eps)
        {
            return std::weak_ordering::less;
        }
        return std::weak_ordering::equivalent;
    }

    constexpr inline bool operator==(Time const& x) const
    {
        return (*this <=> x) == std::weak_ordering::equivalent;
    }

    friend inline std::ostream& operator<<(std::ostream& os, Time const& t)
    {
        auto const precision = os.precision();
        return os << std::setprecision(
                         std::numeric_limits<double>::max_digits10)
                  << t() << std::setprecision(precision);
    }

private:
    MathLib::KahanSum value_;
};

}  // namespace NumLib
