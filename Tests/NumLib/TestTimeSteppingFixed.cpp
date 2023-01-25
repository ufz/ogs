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

TEST(NumLib, TimeSteppingFixed_findDeltatInterval)
{
    double const initial_dt = 1e-12;
    std::vector dts = {initial_dt};
    dts.insert(dts.end(), 10, 1e-1);
    double const t_initial = 0.0;

    {
        std::size_t const expected_interval = 5;
        double time_point = 0.5;
        auto const interval =
            NumLib::findDeltatInterval(t_initial, dts, time_point);
        ASSERT_EQ(expected_interval, interval);
    }
    {
        std::size_t const expected_interval = 5;
        double time_point = 0.5 + std::numeric_limits<double>::epsilon();
        auto const interval =
            NumLib::findDeltatInterval(t_initial, dts, time_point);
        ASSERT_EQ(expected_interval, interval);
    }
    {
        std::size_t const expected_interval = 6;
        double time_point = 0.6 - initial_dt;
        auto const interval_number =
            NumLib::findDeltatInterval(t_initial, dts, time_point);
        INFO("value {} in interval [{}, {})", time_point,
             std::accumulate(begin(dts), begin(dts) + interval_number,
                             initial_dt),
             std::accumulate(begin(dts), begin(dts) + interval_number + 1,
                             initial_dt));
        ASSERT_EQ(expected_interval, interval_number);
    }
    {
        std::size_t const expected_interval =
            std::numeric_limits<std::size_t>::max();
        double time_point = -initial_dt;
        auto const interval_number =
            NumLib::findDeltatInterval(t_initial, dts, time_point);
        ASSERT_EQ(expected_interval, interval_number);
    }
    {
        std::size_t const expected_interval =
            std::numeric_limits<std::size_t>::max();
        double time_point = 1000;
        auto const interval_number =
            NumLib::findDeltatInterval(t_initial, dts, time_point);
        ASSERT_EQ(expected_interval, interval_number);
    }
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

TEST(NumLib, TimeSteppingFixed_incorporateFixedTimesForOutput_2)
{
    double t_initial = 1e-10;
    std::vector<double> dts{};
    dts.insert(dts.end(), 10, 1e-1);
    std::vector<double> fixed_times_for_output{{0.5, 1}};

    auto const expected_time = std::accumulate(dts.begin(), dts.end(), 0.0);

    NumLib::incorporateFixedTimesForOutput(t_initial, expected_time, dts,
                                           fixed_times_for_output);

    // incorporation of time steps doesn't influence the entire simulation time
    EXPECT_NEAR(expected_time, std::accumulate(dts.begin(), dts.end(), 0.0),
                std::numeric_limits<double>::epsilon());

    ASSERT_EQ(1e-1, dts[0]);
    ASSERT_EQ(1e-1, dts[1]);
    ASSERT_EQ(1e-1, dts[2]);
    ASSERT_EQ(1e-1, dts[3]);  // time point 0.4 + 1e-10
    EXPECT_NEAR(1e-1 - t_initial, dts[4],
                std::numeric_limits<double>::epsilon());  // time point 0.5
    EXPECT_EQ(5e-1 + 1e-10 - 0.5, dts[5]);  // time point 0.5 + 1e-10
    EXPECT_EQ(1e-1, dts[6]);                // time point 0.6 + 1e-10
    EXPECT_EQ(1e-1, dts[7]);                // time point 0.7 + 1e-10
    EXPECT_EQ(1e-1, dts[8]);                // time point 0.8 + 1e-10
    EXPECT_EQ(1e-1, dts[9]);                // time point 0.9 + 1e-10
    EXPECT_NEAR(1.0 - (9e-1 + t_initial), dts[10],
                std::numeric_limits<double>::epsilon());  // time point 1.0
    EXPECT_NEAR(
        t_initial, dts[11],
        std::numeric_limits<double>::epsilon());  // time point 1.0 + 1e-10
}

// unit test related to ThermoRichardsMechanics/LiakopoulosHM/liakopoulos.prj
TEST(NumLib, TimeSteppingFixed_incorporateFixedTimesForOutput_3)
{
    double t_initial = 0.0;
    std::vector<double> dts{};

    // <pair> <repeat>10</repeat> <delta_t>1</delta_t> </pair>
    dts.insert(dts.end(), 10, 1);
    // <pair> <repeat>9</repeat> <delta_t>10</delta_t> </pair>
    dts.insert(dts.end(), 9, 10);
    // <pair> <repeat>11</repeat> <delta_t>100</delta_t> </pair>
    dts.insert(dts.end(), 11, 100);
    // <pair> <repeat>3</repeat> <delta_t>200</delta_t> </pair>
    dts.insert(dts.end(), 3, 200);
    // <pair> <repeat>3</repeat> <delta_t>400</delta_t> </pair>
    dts.insert(dts.end(), 3, 400);
    // <pair> <repeat>7</repeat> <delta_t>600</delta_t> </pair>
    dts.insert(dts.end(), 7, 600);

    std::vector<double> fixed_times_for_output{{0.06, 60., 120., 300.0, 600.0,
                                                1200.0, 2400.0, 4800.0, 6000.0,
                                                7200.0}};

    auto const expected_time = std::accumulate(dts.begin(), dts.end(), 0.0);

    NumLib::incorporateFixedTimesForOutput(t_initial, expected_time, dts,
                                           fixed_times_for_output);

    // incorporation of time steps doesn't influence the entire simulation time
    ASSERT_EQ(expected_time, std::accumulate(dts.begin(), dts.end(), 0.0));

    ASSERT_EQ(0.06, dts[0]);
    ASSERT_EQ(1 - 0.06, dts[1]);
    for (std::size_t k = 2; k < 11; ++k)
    {
        ASSERT_EQ(1, dts[k]);
    }
    for (std::size_t k = 11; k < 20; ++k)
    {
        ASSERT_EQ(10, dts[k]);
    }
    ASSERT_EQ(20, dts[20]);  // step size of 100 is split into 20 + 80 because
                             // of output at 120
    ASSERT_EQ(80, dts[21]);
    for (std::size_t k = 22; k < 32; ++k)
    {
        ASSERT_EQ(100, dts[k]);
    }
    for (std::size_t k = 32; k < 35; ++k)
    {
        ASSERT_EQ(200, dts[k]);
    }
    ASSERT_EQ(400, dts[35]);
    ASSERT_EQ(200, dts[36]);  // step size of 400 is split into 200 + 200
                              // because of output at 2400
    ASSERT_EQ(200, dts[37]);
    ASSERT_EQ(400, dts[38]);
    for (std::size_t k = 39; k < 46; ++k)
    {
        ASSERT_EQ(600, dts[k]);
    }
}

// unit test related to
// ogs-ThermoMechanics_CreepBGRa_Verification_m2_3Dload_m2_3Dload
TEST(NumLib, TimeSteppingFixed_incorporateFixedTimesForOutput_4)
{
    double t_initial = 0.0;
    std::vector<double> dts{t_initial};

    // <pair> <repeat>1</repeat> <delta_t>1e-10</delta_t> </pair>
    dts.insert(dts.end(), 1, 1e-10);
    // <pair> <repeat>4</repeat> <delta_t>0.01</delta_t> </pair>
    dts.insert(dts.end(), 4, 1e-2);
    // <pair> <repeat>1</repeat> <delta_t>0.0099999999</delta_t> </pair>
    dts.insert(dts.end(), 1, 0.0099999999);
    // <pair> <repeat>100</repeat> <delta_t>0.01</delta_t> </pair>
    dts.insert(dts.end(), 100, 1e-2);

    std::vector<double> fixed_times_for_output{{0.5, 1.}};

    auto const expected_end_time =
        std::accumulate(dts.begin(), dts.end(), t_initial);

    NumLib::incorporateFixedTimesForOutput(t_initial, expected_end_time, dts,
                                           fixed_times_for_output);

    // incorporation of time steps doesn't influence the entire simulation time
    ASSERT_EQ(expected_end_time,
              std::accumulate(dts.begin(), dts.end(), t_initial));
}

TEST(NumLib, TimeSteppingFixed_incorporateFixedTimesForOutput)
{
    double t_initial = 1.0;
    std::vector<double> dts{{10, 10, 10}};
    std::vector<double> fixed_times_for_output{{9, 12, 28}};

    auto const expected_time = std::accumulate(dts.begin(), dts.end(), 0.0);

    NumLib::incorporateFixedTimesForOutput(t_initial, expected_time, dts,
                                           fixed_times_for_output);

    ASSERT_EQ(expected_time, std::accumulate(dts.begin(), dts.end(), 0.0));
    ASSERT_EQ(8.0, dts[0]);
    ASSERT_EQ(2.0, dts[1]);
    ASSERT_EQ(1.0, dts[2]);
    ASSERT_EQ(9.0, dts[3]);
    ASSERT_EQ(7.0, dts[4]);
    ASSERT_EQ(3.0, dts[5]);
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
