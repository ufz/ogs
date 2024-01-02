/**
 * \file
 * \author Norihiro Watanabe
 * \date   2012-08-03
 *
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
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

class NumLibTimeSteppingFixed_FindDeltaT : public ::testing::Test
{
    std::vector<double> dts;

protected:
    double const initial_dt = 1e-12;

    void SetUp() override
    {
        dts = {initial_dt};
        dts.insert(dts.end(), 10, 1e-1);
    }

    std::size_t find(double const t)
    {
        double const t_initial = 0.0;
        return NumLib::findDeltatInterval(t_initial, dts, t);
    }
};

TEST_F(NumLibTimeSteppingFixed_FindDeltaT, TestPointWithinMiddleInterval)
{
    EXPECT_EQ(5, find(0.5));
}

TEST_F(NumLibTimeSteppingFixed_FindDeltaT, TestPointAtBoundary)
{
    EXPECT_EQ(5, find(0.5 + std::numeric_limits<double>::epsilon()));
}

TEST_F(NumLibTimeSteppingFixed_FindDeltaT, TestPointJustBeforeNextInterval)
{
    EXPECT_EQ(6, find(0.6 - initial_dt));
}

TEST_F(NumLibTimeSteppingFixed_FindDeltaT, TestPointBeforeInitialTime)
{
    EXPECT_EQ(std::numeric_limits<std::size_t>::max(), find(-initial_dt));
}

TEST_F(NumLibTimeSteppingFixed_FindDeltaT, TestPointBeyondEndTime)
{
    EXPECT_EQ(std::numeric_limits<std::size_t>::max(), find(1000));
}

class NumLibTimeSteppingFixed_TimeSteps : public ::testing::Test
{
    std::vector<int> const dummy_number_iterations = {0, 0, 0, 0};

protected:
    void checkExpectedTimes(NumLib::FixedTimeStepping& algorithm,
                            std::vector<double> const expected_times) const
    {
        std::vector<double> times =
            timeStepping(algorithm, dummy_number_iterations, {});

        ASSERT_EQ(expected_times.size(), times.size());
        ASSERT_ARRAY_NEAR(expected_times, times, expected_times.size(),
                          std::numeric_limits<double>::epsilon());
    }
};

TEST_F(NumLibTimeSteppingFixed_TimeSteps, HomogeneousDt)
{
    NumLib::FixedTimeStepping algorithm(1, 31, 10);

    checkExpectedTimes(algorithm, {1, 11, 21, 31});
}

TEST_F(NumLibTimeSteppingFixed_TimeSteps, DtVectorEqualSum)
{
    std::vector<std::pair<std::size_t, double>> const repeat_dt = {{3, 10}};
    NumLib::FixedTimeStepping algorithm(1, 31, repeat_dt, {});

    checkExpectedTimes(algorithm, {1, 11, 21, 31});
}

TEST_F(NumLibTimeSteppingFixed_TimeSteps, DtVectorLessThanSum)
{
    std::vector<std::pair<std::size_t, double>> const repeat_dt = {
        {1, 5}, {1, 10}, {1, 20}};
    NumLib::FixedTimeStepping algorithm(1, 31, repeat_dt, {});

    checkExpectedTimes(algorithm, {1, 6, 16, 31});
}

TEST_F(NumLibTimeSteppingFixed_TimeSteps, DtVectorGreaterThanSum)
{
    std::vector<std::pair<std::size_t, double>> const repeat_dt = {{1, 5},
                                                                   {3, 10}};
    NumLib::FixedTimeStepping algorithm(1, 31, repeat_dt, {});

    checkExpectedTimes(algorithm, {1, 6, 16, 26, 31});
}

TEST(NumLibTimeSteppingFixed_FixedOutputTimes, incorporateFixedTimesForOutput_2)
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
TEST(NumLibTimeSteppingFixed_FixedOutputTimes, incorporateFixedTimesForOutput_3)
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
TEST(NumLibTimeSteppingFixed_FixedOutputTimes, incorporateFixedTimesForOutput_4)
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

TEST(NumLibTimeSteppingFixed_FixedOutputTimes, incorporateFixedTimesForOutput)
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

TEST(NumLibTimeSteppingFixed_FixedOutputTimes,
     incorporateFixedTimesForOutput_Matching)
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

TEST(NumLibTimeSteppingFixed_FixedOutputTimes,
     OutputTimeBeforeSimulationStartTime)
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

TEST(NumLibTimeSteppingFixed_FixedOutputTimes, OutputTimeAfterSimulationEndTime)
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
