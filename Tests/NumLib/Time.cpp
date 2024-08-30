/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "NumLib/TimeStepping/Time.h"

#include <gtest/gtest.h>

#include <numeric>
#include <vector>

TEST(NumLibTime, AddSameIncrementToEquivalentTimes)
{
    // two time points 'set' (set to one) and 'numeric' computed by summation
    // such that the same time point as 'set' is reached (but the correction may
    // be different)
    NumLib::Time time_numeric(0.0);
    NumLib::Time time_set(1.0);
    double const delta_t = 1e-8;
    for (int i = 0; i < 1e8; ++i)
    {
        time_numeric += delta_t;
    }
    EXPECT_EQ(time_set, time_numeric);
}

TEST(NumLibTime, AddSameIncrementsInDifferentOrder)
{
    NumLib::Time time_a(0.0);
    NumLib::Time time_b(0.0);

    time_a += 1.0;
    time_a += 1e17;
    time_a -= 1.0;
    time_a -= 1e17;

    time_b += 1e17;
    time_b -= 1e17;
    time_b -= 1.0;
    time_b += 1.0;
    EXPECT_EQ(time_a, time_b);
}

TEST(NumLibTime, SmallIncrements)
{
    double const delta_t = 1e-8;
    NumLib::Time current_time(1.0);
    for (int i = 1; i < 1e8 + 1; ++i)
    {
        current_time += delta_t;
        ASSERT_EQ(NumLib::Time{1.0 + i * delta_t}, current_time);
    }
    ASSERT_EQ(NumLib::Time{2.0}, current_time);
}

TEST(NumLibTime, SmallAndLargeIncrements)
{
    double const delta_t_small = 1e-8;
    double const delta_t_large = 10.3;
    NumLib::Time current_time(1.0);
    for (int i = 1; i < 1e7 + 1; ++i)
    {
        current_time += delta_t_small;
        ASSERT_EQ(NumLib::Time{1.0 + i * delta_t_small}, current_time);
    }
    ASSERT_EQ(NumLib::Time{1.1}, current_time);
    for (int i = 1; i < 10; ++i)
    {
        current_time += delta_t_large;
        ASSERT_EQ(NumLib::Time{1.1 + i * delta_t_large}, current_time);
    }
    ASSERT_EQ(NumLib::Time{93.8}, current_time);
    for (int i = 1; i < 1e7 + 1; ++i)
    {
        current_time += delta_t_small;
        ASSERT_EQ(NumLib::Time{93.8 + i * delta_t_small}, current_time);
    }
    ASSERT_EQ(NumLib::Time{93.9}, current_time);
}

TEST(NumLibTime, TH2M_TH_idealGasLawBenchmark)
{
    double const delta_t = 1e-1;
    NumLib::Time current_time(0.0);
    for (int i = 0; i < 100; ++i)
    {
        current_time = current_time + delta_t;
        ASSERT_EQ(NumLib::Time{(i + 1) * delta_t}, current_time);
    }
    ASSERT_EQ(NumLib::Time{10.0}, current_time);
}
