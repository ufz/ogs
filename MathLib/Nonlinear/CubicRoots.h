/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */
#pragma once

#include <algorithm>
#include <limits>
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

    // Method to solve the cubic equation using Boost
    std::vector<double> solve()
    {
        // Solve using Boost's cubic_roots
        auto roots = boost::math::tools::cubic_roots<double>(a_, b_, c_, d_);

        // Sort array
        std::sort(roots.begin(), roots.end());

        std::vector<double> roots_vector;
        for (const auto& root : roots)
        {
            if (!std::isnan(root))
            {
                roots_vector.push_back(root);
            }
        }
        return roots_vector;
    }

    // Method that returns the smallest positive real root
    double smallestPositiveRealRoot()
    {
        auto roots = solve();
        double min_positive_root = std::numeric_limits<double>::infinity();

        // Find the smallest positive real root
        for (const auto& root : roots)
        {
            if (root > 0 && root < min_positive_root)
            {
                min_positive_root = root;
            }
        }

        // If no positive root was found, return NaN
        if (min_positive_root == std::numeric_limits<double>::infinity())
        {
            return std::numeric_limits<double>::quiet_NaN();
        }

        return min_positive_root;
    }

private:
    double a_, b_, c_, d_;  // Coefficients of the cubic equation
};

}  // namespace MathLib
