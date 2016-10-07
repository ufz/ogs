/**
 * @file TestPiecewiseLinearInterpolation.cpp
 * @author Thomas Fischer
 * @date Feb 12, 2013
 * @brief
 *
 * @copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

// stl
#include <limits>

// google test
#include "gtest/gtest.h"

#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"

TEST(MathLibInterpolationAlgorithms, PiecewiseLinearInterpolation)
{
    const std::size_t size(1000);
    std::vector<double> supp_pnts, values;
    for (std::size_t k(0); k < size; ++k)
    {
        supp_pnts.push_back(static_cast<double>(k));
        if (k % 2 == 0)
        {
            values.push_back(0.0);
        }
        else
        {
            values.push_back(1.0);
        }
    }

    MathLib::PiecewiseLinearInterpolation interpolation{std::move(supp_pnts),
                                                        std::move(values)};
    // Interpolation
    for (std::size_t k(0); k < size - 1; ++k)
    {
        ASSERT_NEAR(0.5, interpolation.getValue(k + 0.5),
                    std::numeric_limits<double>::epsilon());
    }

    for (std::size_t k(0); k < size - 1; ++k)
    {
        if (k % 2 == 0)
        {
            ASSERT_NEAR(0.25, interpolation.getValue(k + 0.25),
                        std::numeric_limits<double>::epsilon());
            ASSERT_NEAR(0.75, interpolation.getValue(k + 0.75),
                        std::numeric_limits<double>::epsilon());
        }
        else
        {
            ASSERT_NEAR(0.75, interpolation.getValue(k + 0.25),
                        std::numeric_limits<double>::epsilon());
            ASSERT_NEAR(0.25, interpolation.getValue(k + 0.75),
                        std::numeric_limits<double>::epsilon());
        }
    }

    // Extrapolation
    ASSERT_NEAR(0.0, interpolation.getValue(-0.5),
                std::numeric_limits<double>::epsilon());
    // Extrapolation
    ASSERT_NEAR(1.0, interpolation.getValue(size - 0.5),
                std::numeric_limits<double>::epsilon());
}

TEST(MathLibInterpolationAlgorithms,
     PiecewiseLinearInterpolationSupportPntsInReverseOrder)
{
    const std::size_t size(1000);
    std::vector<double> supp_pnts, values;
    for (std::size_t k(0); k < size; ++k)
    {
        supp_pnts.push_back(static_cast<double>(size - 1 - k));
        if (k % 2 == 0)
        {
            values.push_back(1.0);
        }
        else
        {
            values.push_back(0.0);
        }
    }

    MathLib::PiecewiseLinearInterpolation interpolation{std::move(supp_pnts),
                                                        std::move(values)};
    // Interpolation
    for (std::size_t k(0); k < size - 1; ++k)
    {
        ASSERT_NEAR(0.5, interpolation.getValue(k + 0.5),
                    std::numeric_limits<double>::epsilon());
    }

    for (std::size_t k(0); k < size - 1; ++k)
    {
        if (k % 2 == 0)
        {
            ASSERT_NEAR(0.25, interpolation.getValue(k + 0.25),
                        std::numeric_limits<double>::epsilon());
            ASSERT_NEAR(0.75, interpolation.getValue(k + 0.75),
                        std::numeric_limits<double>::epsilon());
        }
        else
        {
            ASSERT_NEAR(0.75, interpolation.getValue(k + 0.25),
                        std::numeric_limits<double>::epsilon());
            ASSERT_NEAR(0.25, interpolation.getValue(k + 0.75),
                        std::numeric_limits<double>::epsilon());
        }
    }

    // Extrapolation
    ASSERT_NEAR(0.0, interpolation.getValue(-0.5),
                std::numeric_limits<double>::epsilon());
    // Extrapolation
    ASSERT_NEAR(1.0, interpolation.getValue(size - 0.5),
                std::numeric_limits<double>::epsilon());
}

TEST(MathLibInterpolationAlgorithms, PiecewiseLinearInterpolationDerivative)
{
    const std::size_t size(1000);
    std::vector<double> supp_pnts, values;
    for (std::size_t k(0); k < size; ++k)
    {
        supp_pnts.push_back(static_cast<double>(k));
        values.push_back(k * k);
    }

    MathLib::PiecewiseLinearInterpolation interpolation{std::move(supp_pnts),
                                                        std::move(values)};
    // Interpolation
    for (std::size_t k(0); k < size - 1; ++k)
    {
        ASSERT_NEAR(1 + 2 * k, interpolation.GetDerivative(k + 0.5),
                    std::numeric_limits<double>::epsilon());
    }

    // Extrapolation
    ASSERT_NEAR(0, interpolation.GetDerivative(-1),
                std::numeric_limits<double>::epsilon());
    // Extrapolation
    ASSERT_NEAR(0, interpolation.GetDerivative(1001),
                std::numeric_limits<double>::epsilon());
}
