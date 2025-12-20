// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <spdlog/fmt/ostr.h>

#include <iomanip>
#include <ostream>
#include <range/v3/range/concepts.hpp>
#include <range/v3/range/traits.hpp>

namespace MathLib
{
class KahanSum;
}
template <>
struct fmt::formatter<MathLib::KahanSum> : fmt::ostream_formatter
{
};

namespace MathLib
{

class KahanSum
{
public:
    explicit constexpr KahanSum(double const value = 0) : value_(value) {}

    explicit constexpr KahanSum(ranges::range auto const& range) : value_(0)
    {
        for (auto const v : range)
        {
            *this += v;
        }
    }

    constexpr KahanSum operator+(double const increment) const
    {
        KahanSum result = *this;
        return result += increment;
    }

    constexpr KahanSum operator-(double const increment) const
    {
        KahanSum result = *this;
        return result += -increment;
    }

    constexpr KahanSum& operator-=(double const increment)
    {
        return *this += -increment;
    }

    constexpr KahanSum& operator+=(double const increment)
    {
        double const y = increment - correction_;
        double const t = value_ + y;
        correction_ = (t - value_) - y;
        value_ = t;
        return *this;
    }

    constexpr double value() const { return value_; }
    constexpr double operator()() const { return value_; }

    friend inline std::ostream& operator<<(std::ostream& os, KahanSum const& x)
    {
        auto const precision = os.precision();
        return os << std::setprecision(
                         std::numeric_limits<double>::max_digits10)
                  << x.value() << " (± " << x.correction_ << ')'
                  << std::setprecision(precision);
    }

private:
    double value_;
    double correction_ = 0.;
};

}  // namespace MathLib
