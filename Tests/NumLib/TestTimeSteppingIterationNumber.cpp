/**
 * \author Norihiro Watanabe
 * \date   2012-08-03
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <gtest/gtest.h>

#include <utility>
#include <vector>

#include "BaseLib/Logging.h"

#include "NumLib/TimeStepping/TimeStep.h"
#include "NumLib/TimeStepping/Algorithms/IterationNumberBasedTimeStepping.h"

#include "Tests/TestTools.h"
#include "TimeSteppingTestingTools.h"

TEST(NumLib, TimeSteppingIterationNumberBased1)
{
    std::vector<int> iter_times_vector = {0, 3, 5, 7};
    std::vector<double> multiplier_vector = {2.0, 1.0, 0.5, 0.25};
    NumLib::IterationNumberBasedTimeStepping alg(
        1, 31, 1, 10, 1, std::move(iter_times_vector),
        std::move(multiplier_vector), {});

    const double solution_error = 0.;

    ASSERT_TRUE(alg.next(solution_error, 1));
    NumLib::TimeStep ts = alg.getTimeStep();
    ASSERT_EQ(1u, ts.steps());
    ASSERT_EQ(1., ts.previous());
    ASSERT_EQ(2., ts.current());
    ASSERT_EQ(1., ts.dt());
    ASSERT_TRUE(alg.accepted());

    ASSERT_TRUE(alg.next(solution_error, 1));

    ASSERT_TRUE(alg.next(solution_error, 3));
    ts = alg.getTimeStep();
    ASSERT_EQ(3u, ts.steps());
    ASSERT_EQ(4., ts.previous());
    ASSERT_EQ(6., ts.current());
    ASSERT_EQ(2., ts.dt());
    ASSERT_TRUE(alg.accepted());

    ASSERT_TRUE(alg.next(solution_error, 5));
    ts = alg.getTimeStep();
    ASSERT_EQ(4u, ts.steps());
    ASSERT_EQ(6., ts.previous());
    ASSERT_EQ(7., ts.current());
    ASSERT_EQ(1., ts.dt());
    ASSERT_TRUE(alg.accepted());

    ASSERT_TRUE(alg.next(solution_error, 7));
    ts = alg.getTimeStep();
    ASSERT_EQ(5u, ts.steps());
    ASSERT_EQ(7., ts.previous());
    ASSERT_EQ(8., ts.current());
    ASSERT_EQ(1., ts.dt());
    ASSERT_TRUE(alg.accepted());

    ASSERT_TRUE(alg.next(solution_error, 8 /* exceed maximum */));
    ts = alg.getTimeStep();
    ASSERT_EQ(6u, ts.steps());
    ASSERT_EQ(8., ts.previous());
    ASSERT_EQ(9, ts.current());
    ASSERT_EQ(1., ts.dt());
    ASSERT_TRUE(alg.accepted());

    ASSERT_TRUE(alg.next(solution_error, 4));
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
    NumLib::IterationNumberBasedTimeStepping alg(
        1, 31, 1, 10, 1, std::move(iter_times_vector),
        std::move(multiplier_vector), {});

    std::vector<int> nr_iterations = {0, 2, 2, 2, 4, 6, 8, 4, 1};
    const std::vector<double> expected_vec_t = {1,  2,  4,  8,  16,
                                                24, 28, 29, 30, 31};

    std::vector<double> vec_t = timeStepping(alg, nr_iterations);

    ASSERT_EQ(expected_vec_t.size(), vec_t.size());
    ASSERT_EQ(0u, alg.getNumberOfRepeatedSteps());
    ASSERT_ARRAY_NEAR(expected_vec_t, vec_t, expected_vec_t.size(), std::numeric_limits<double>::epsilon());
}

TEST(NumLib, TimeSteppingIterationNumberBased2FixedOutputTimes)
{
    std::vector<int> iter_times_vector = {0, 3, 5, 7};
    std::vector<double> multiplier_vector = {2.0, 1.0, 0.5, 0.25};
    std::vector<double> fixed_output_times = {5, 20};
    NumLib::IterationNumberBasedTimeStepping alg(
        1, 31, 1, 10, 1, std::move(iter_times_vector),
        std::move(multiplier_vector), std::move(fixed_output_times));

    std::vector<int> nr_iterations = {0, 2, 2, 2, 4, 6, 8, 4, 1, 1, 1, 1, 1};
    const std::vector<double> expected_vec_t = {1,  2,  4,  5,  7,  9,  10,
                                                11, 12, 14, 18, 20, 24, 31};

    std::vector<double> vec_t = timeStepping(alg, nr_iterations);

    EXPECT_EQ(expected_vec_t.size(), vec_t.size());
    ASSERT_EQ(0u, alg.getNumberOfRepeatedSteps());
    ASSERT_ARRAY_NEAR(expected_vec_t, vec_t, expected_vec_t.size(),
                      std::numeric_limits<double>::epsilon());
}
