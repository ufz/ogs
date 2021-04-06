/**
 * \author Norihiro Watanabe
 * \date   2012-08-03
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <gtest/gtest.h>

#include <tuple>
#include <utility>
#include <vector>

#include "NumLib/TimeStepping/Algorithms/IterationNumberBasedTimeStepping.h"
#include "NumLib/TimeStepping/TimeStep.h"
#include "Tests/TestTools.h"
#include "TimeSteppingTestingTools.h"

TEST(NumLib, TimeSteppingIterationNumberBased1)
{
    std::vector<int> iter_times_vector = {0, 3, 5, 7};
    std::vector<double> multiplier_vector = {2.0, 1.0, 0.5, 0.25};
    NumLib::IterationNumberBasedTimeStepping alg(1, 31, 1, 10, 1,
                                                 std::move(iter_times_vector),
                                                 std::move(multiplier_vector));

    const double solution_error = 0.;
    const double end_time = alg.end();
    auto [step_accepted, timestepper_dt] = alg.next(solution_error, 1);
    ASSERT_TRUE(step_accepted);
    timestepper_dt = (alg.getTimeStep().current() + timestepper_dt > end_time)
                         ? end_time - alg.getTimeStep().current()
                         : timestepper_dt;
    alg.resetCurrentTimeStep(timestepper_dt);
    NumLib::TimeStep ts = alg.getTimeStep();
    ASSERT_EQ(1u, ts.steps());
    ASSERT_EQ(1., ts.previous());
    ASSERT_EQ(2., ts.current());
    ASSERT_EQ(1., ts.dt());
    ASSERT_TRUE(alg.accepted());

    auto [step_accepted1, timestepper_dt1] = alg.next(solution_error, 1);
    ASSERT_TRUE(step_accepted1);
    timestepper_dt1 = (alg.getTimeStep().current() + timestepper_dt1 > end_time)
                          ? end_time - alg.getTimeStep().current()
                          : timestepper_dt1;
    alg.resetCurrentTimeStep(timestepper_dt1);

    auto [step_accepted2, timestepper_dt2] = alg.next(solution_error, 3);
    ASSERT_TRUE(step_accepted2);
    timestepper_dt2 = (alg.getTimeStep().current() + timestepper_dt2 > end_time)
                          ? end_time - alg.getTimeStep().current()
                          : timestepper_dt2;
    alg.resetCurrentTimeStep(timestepper_dt2);
    ts = alg.getTimeStep();
    ASSERT_EQ(3u, ts.steps());
    ASSERT_EQ(4., ts.previous());
    ASSERT_EQ(6., ts.current());
    ASSERT_EQ(2., ts.dt());
    ASSERT_TRUE(alg.accepted());

    auto [step_accepted3, timestepper_dt3] = alg.next(solution_error, 5);
    ASSERT_TRUE(step_accepted3);
    timestepper_dt3 = (alg.getTimeStep().current() + timestepper_dt3 > end_time)
                          ? end_time - alg.getTimeStep().current()
                          : timestepper_dt3;
    alg.resetCurrentTimeStep(timestepper_dt3);
    ts = alg.getTimeStep();
    ASSERT_EQ(4u, ts.steps());
    ASSERT_EQ(6., ts.previous());
    ASSERT_EQ(7., ts.current());
    ASSERT_EQ(1., ts.dt());
    ASSERT_TRUE(alg.accepted());

    auto [step_accepted4, timestepper_dt4] = alg.next(solution_error, 7);
    ASSERT_TRUE(step_accepted4);
    timestepper_dt4 = (alg.getTimeStep().current() + timestepper_dt4 > end_time)
                          ? end_time - alg.getTimeStep().current()
                          : timestepper_dt4;
    alg.resetCurrentTimeStep(timestepper_dt4);
    ts = alg.getTimeStep();
    ASSERT_EQ(5u, ts.steps());
    ASSERT_EQ(7., ts.previous());
    ASSERT_EQ(8., ts.current());
    ASSERT_EQ(1., ts.dt());
    ASSERT_TRUE(alg.accepted());

    auto [step_accepted5, timestepper_dt5] = alg.next(solution_error, 8);
    ASSERT_TRUE(step_accepted5);
    timestepper_dt5 = (alg.getTimeStep().current() + timestepper_dt5 > end_time)
                          ? end_time - alg.getTimeStep().current()
                          : timestepper_dt5;
    alg.resetCurrentTimeStep(timestepper_dt5);
    ts = alg.getTimeStep();
    ASSERT_EQ(6u, ts.steps());
    ASSERT_EQ(8., ts.previous());
    ASSERT_EQ(9, ts.current());
    ASSERT_EQ(1., ts.dt());
    ASSERT_TRUE(alg.accepted());

    auto [step_accepted6, timestepper_dt6] = alg.next(solution_error, 4);
    ASSERT_TRUE(step_accepted6);
    timestepper_dt6 = (alg.getTimeStep().current() + timestepper_dt6 > end_time)
                          ? end_time - alg.getTimeStep().current()
                          : timestepper_dt6;
    alg.resetCurrentTimeStep(timestepper_dt6);
    ts = alg.getTimeStep();
    ASSERT_EQ(7u, ts.steps());
    ASSERT_EQ(9., ts.previous());
    ASSERT_EQ(10, ts.current());
    ASSERT_EQ(1., ts.dt());
    ASSERT_TRUE(alg.accepted());
}

TEST(NumLib, TimeSteppingIterationNumberBased2)
{
    std::vector<int> iter_times_vector = {0, 3, 5, 7};
    std::vector<double> multiplier_vector = {2.0, 1.0, 0.5, 0.25};
    NumLib::IterationNumberBasedTimeStepping alg(1, 31, 1, 10, 1,
                                                 std::move(iter_times_vector),
                                                 std::move(multiplier_vector));

    std::vector<int> nr_iterations = {0, 2, 2, 2, 4, 6, 8, 4, 1};
    const std::vector<double> expected_vec_t = {1,  2,  4,  8,  16,
                                                24, 28, 29, 30, 31};

    std::vector<double> vec_t = timeStepping(alg, nr_iterations, {});

    ASSERT_EQ(expected_vec_t.size(), vec_t.size());
    ASSERT_ARRAY_NEAR(expected_vec_t, vec_t, expected_vec_t.size(),
                      std::numeric_limits<double>::epsilon());
}

TEST(NumLib, TimeSteppingIterationNumberBased2FixedOutputTimes)
{
    std::vector<int> iter_times_vector = {0, 3, 5, 7};
    std::vector<double> multiplier_vector = {2.0, 1.0, 0.5, 0.25};
    std::vector<double> fixed_output_times = {5, 20};
    NumLib::IterationNumberBasedTimeStepping alg(1, 31, 1, 10, 1,
                                                 std::move(iter_times_vector),
                                                 std::move(multiplier_vector));

    std::vector<int> nr_iterations = {0, 2, 2, 2, 4, 6, 8, 4, 1, 1, 1, 1, 1};
    const std::vector<double> expected_vec_t = {1,  2,  4,  5,  7,  9,  10,
                                                11, 12, 14, 18, 20, 24, 31};

    std::vector<double> vec_t =
        timeStepping(alg, nr_iterations, fixed_output_times);

    EXPECT_EQ(expected_vec_t.size(), vec_t.size());
    ASSERT_ARRAY_NEAR(expected_vec_t, vec_t, expected_vec_t.size(),
                      std::numeric_limits<double>::epsilon());
}
