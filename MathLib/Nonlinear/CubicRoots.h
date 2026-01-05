// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <algorithm>
#include <limits>
#include <range/v3/algorithm/min.hpp>
#include <range/v3/algorithm/sort.hpp>
#include <range/v3/range/conversion.hpp>
#include <range/v3/view/filter.hpp>
#include <vector>

#include "BaseLib/Error.h"
#include "cubic_roots.hpp"
namespace MathLib
{

class CubicSolver
{
public:
    CubicSolver(const double a, const double b, const double c, const double d)
        : a_(a), b_(b), c_(c), d_(d)
    {
        if (std::abs(a_) < 1e-9)
        {
            OGS_FATAL("'a' must be non-zero for a cubic equation.");
        }
    }

    std::vector<double> solve()
    {
        std::array<double, 3> const roots =
            boost::math::tools::cubic_roots<double>(a_, b_, c_, d_);

        std::vector<double> adjusted_roots;
        adjusted_roots.reserve(3);

        double last_valid = std::numeric_limits<double>::quiet_NaN();

        for (auto root : roots)
        {
            if (std::isnan(root))
            {
                if (!std::isnan(last_valid))
                {
                    adjusted_roots.push_back(last_valid);
                }
                // If we get NaN before any valid root, we just skip it
            }
            else
            {
                adjusted_roots.push_back(root);
                last_valid = root;
            }
        }

        ranges::sort(adjusted_roots);
        return adjusted_roots;
    }
    // Method that returns the smallest positive real root
    double smallestPositiveRealRoot()
    {
        std::vector<double> const roots = solve();

        auto positive_roots =
            roots | ranges::views::filter([](double root) { return root > 0; });

        // If no positive root exists, return NaN
        if (ranges::empty(positive_roots))
        {
            return std::numeric_limits<double>::quiet_NaN();
        }

        return ranges::min(positive_roots);
    }

private:
    double a_, b_, c_, d_;  // Coefficients of the cubic equation
};

}  // namespace MathLib
