/**
 * \file
 * \author Norihiro Watanabe
 * \date   2012-08-03
 *
 * \copyright
 * Copyright (c) 2012-2023, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <gtest/gtest.h>

#include <numeric>
#include <vector>

#include "NumLib/TimeStepping/Algorithms/CreateFixedTimeStepping.h"
#include "NumLib/TimeStepping/Algorithms/FixedTimeStepping.h"
#include "NumLib/TimeStepping/TimeStep.h"
#include "Tests/TestTools.h"
#include "TimeSteppingTestingTools.h"

namespace NumLib
{
extern void incorporateFixedTimesForOutput(
    double t_initial, double t_end, std::vector<double>& timesteps,
    std::vector<double> const& fixed_times_for_output);
}

TEST(NumLib, TimeSteppingFixed)
{
    std::vector<int> const dummy_number_iterations = {0, 0, 0, 0};
    // homogeneous dt
    {
        NumLib::FixedTimeStepping fixed(1, 31, 10);
        const std::vector<double> expected_vec_t = {1, 11, 21, 31};

        std::vector<double> vec_t =
            timeStepping(fixed, dummy_number_iterations, {});

        ASSERT_EQ(expected_vec_t.size(), vec_t.size());
        ASSERT_ARRAY_NEAR(expected_vec_t, vec_t, expected_vec_t.size(),
                          std::numeric_limits<double>::epsilon());
    }

    // dt vector (t_end == t0 + sum(dt))
    {
        const std::vector<double> fixed_dt = {10, 10, 10};
        NumLib::FixedTimeStepping fixed(1, 31, fixed_dt);
        const std::vector<double> expected_vec_t = {1, 11, 21, 31};

        std::vector<double> vec_t =
            timeStepping(fixed, dummy_number_iterations, {});

        ASSERT_EQ(expected_vec_t.size(), vec_t.size());
        ASSERT_ARRAY_NEAR(expected_vec_t, vec_t, expected_vec_t.size(),
                          std::numeric_limits<double>::epsilon());
    }

    // dt vector (t_end < t0 + sum(dt))
    {
        const std::vector<double> fixed_dt = {5, 10, 20};
        NumLib::FixedTimeStepping fixed(1, 31, fixed_dt);
        const std::vector<double> expected_vec_t = {1, 6, 16, 31};

        std::vector<double> vec_t =
            timeStepping(fixed, dummy_number_iterations, {});

        ASSERT_EQ(expected_vec_t.size(), vec_t.size());
        ASSERT_ARRAY_NEAR(expected_vec_t, vec_t, expected_vec_t.size(),
                          std::numeric_limits<double>::epsilon());
    }

    // dt vector (t_end > t0 + sum(dt))
    {
        const std::vector<double> fixed_dt = {5, 10, 10};
        NumLib::FixedTimeStepping fixed(1, 31, fixed_dt);
        const std::vector<double> expected_vec_t = {1, 6, 16, 26};

        std::vector<double> vec_t =
            timeStepping(fixed, dummy_number_iterations, {});

        ASSERT_EQ(expected_vec_t.size(), vec_t.size());
        ASSERT_ARRAY_NEAR(expected_vec_t, vec_t, expected_vec_t.size(),
                          std::numeric_limits<double>::epsilon());
    }
}

TEST(NumLib, TimeSteppingFixed_incorporateFixedTimesForOutput)
{
    double t_initial = 1.0;
    std::vector<double> timesteps{{10, 10, 10}};
    std::vector<double> fixed_times_for_output{{9, 12, 28}};

    auto const expected_time =
        std::accumulate(timesteps.begin(), timesteps.end(), 0.0);

    NumLib::incorporateFixedTimesForOutput(t_initial, expected_time, dts,
                                           fixed_times_for_output);

    ASSERT_EQ(expected_time,
              std::accumulate(timesteps.begin(), timesteps.end(), 0.0));
    ASSERT_EQ(8.0, timesteps[0]);
    ASSERT_EQ(2.0, timesteps[1]);
    ASSERT_EQ(1.0, timesteps[2]);
    ASSERT_EQ(9.0, timesteps[3]);
    ASSERT_EQ(7.0, timesteps[4]);
    ASSERT_EQ(3.0, timesteps[5]);
}

TEST(NumLib, TimeSteppingFixed_incorporateFixedTimesForOutput_Matching)
{
    double t_initial = 1.0;
    std::vector<double> timesteps{{10, 10, 10}};
    std::vector<double> fixed_times_for_output{{9, 11, 31}};

    auto const expected_time =
        std::accumulate(timesteps.begin(), timesteps.end(), 0.0);

    NumLib::incorporateFixedTimesForOutput(t_initial, expected_time, timesteps,
                                           fixed_times_for_output);

    ASSERT_EQ(expected_time,
              std::accumulate(timesteps.begin(), timesteps.end(), 0.0));
    ASSERT_EQ(8.0, timesteps[0]);
    ASSERT_EQ(2.0, timesteps[1]);
    ASSERT_EQ(10.0, timesteps[2]);
    ASSERT_EQ(10.0, timesteps[3]);
}

TEST(
    NumLib,
    TimeSteppingFixed_incorporateFixedTimesForOutput_OutputTimeBeforeSimulationStartTime)
{
    double t_initial = 10.0;
    std::vector<double> timesteps{{10, 10, 10}};
    std::vector<double> fixed_times_for_output{{9, 12, 28}};

    auto const expected_time =
        std::accumulate(timesteps.begin(), timesteps.end(), t_initial);

    NumLib::incorporateFixedTimesForOutput(t_initial, expected_time, timesteps,
                                           fixed_times_for_output);
    ASSERT_EQ(expected_time,
              std::accumulate(timesteps.begin(), timesteps.end(), t_initial));
    ASSERT_EQ(2.0, timesteps[0]);
    ASSERT_EQ(8.0, timesteps[1]);
    ASSERT_EQ(8.0, timesteps[2]);
    ASSERT_EQ(2.0, timesteps[3]);
    ASSERT_EQ(10.0, timesteps[4]);
}

TEST(
    NumLib,
    TimeSteppingFixed_incorporateFixedTimesForOutput_OutputTimeAfterSimulationEndTime)
{
    double t_initial = 1.0;
    std::vector<double> timesteps{{10, 10, 10}};
    std::vector<double> fixed_times_for_output{{9, 12, 28, 33}};

    auto const expected_time =
        std::accumulate(timesteps.begin(), timesteps.end(), t_initial);

    NumLib::incorporateFixedTimesForOutput(t_initial, expected_time, timesteps,
                                           fixed_times_for_output);
    ASSERT_EQ(expected_time,
              std::accumulate(timesteps.begin(), timesteps.end(), t_initial));
    ASSERT_EQ(8.0, timesteps[0]);
    ASSERT_EQ(2.0, timesteps[1]);
    ASSERT_EQ(1.0, timesteps[2]);
    ASSERT_EQ(9.0, timesteps[3]);
    ASSERT_EQ(7.0, timesteps[4]);
    ASSERT_EQ(3.0, timesteps[5]);
}
