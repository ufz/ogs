/**
 * \author Norihiro Watanabe
 * \date   2012-08-03
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <gtest/gtest.h>

#include <utility>
#include <vector>

#include <logog/include/logog.hpp>

#include "BaseLib/DebugTools.h"
#include "NumLib/TimeStepping/TimeStep.h"
#include "NumLib/TimeStepping/Algorithms/IterationNumberBasedAdaptiveTimeStepping.h"

#include "Tests/TestTools.h"
#include "TimeSteppingTestingTools.h"

TEST(NumLib, TimeSteppingIterationNumberBased1)
{
    std::vector<std::size_t> iter_times_vector = {0, 3, 5, 7};
    std::vector<double> multiplier_vector = {2.0, 1.0, 0.5, 0.25};
    NumLib::IterationNumberBasedAdaptiveTimeStepping alg(1, 31, 1, 10, 1, iter_times_vector, multiplier_vector);

    const double solution_error = 0.;

    ASSERT_TRUE(alg.next(solution_error));  // t=2, dt=1
    NumLib::TimeStep ts = alg.getTimeStep();
    ASSERT_EQ(1u, ts.steps());
    ASSERT_EQ(1., ts.previous());
    ASSERT_EQ(2., ts.current());
    ASSERT_EQ(1., ts.dt());
    ASSERT_TRUE(alg.accepted());

    ASSERT_TRUE(alg.next(solution_error));  // t=4, dt=2

    // dt*=2
    alg.setIterationNumber(3);
    ASSERT_TRUE(alg.next(solution_error));  // t=8, dt=4
    ts = alg.getTimeStep();
    ASSERT_EQ(3u, ts.steps());
    ASSERT_EQ(4., ts.previous());
    ASSERT_EQ(8., ts.current());
    ASSERT_EQ(4., ts.dt());
    ASSERT_TRUE(alg.accepted());

    // dt*=1
    alg.setIterationNumber(5);
    ASSERT_TRUE(alg.next(solution_error));  // t=12, dt=4
    ts = alg.getTimeStep();
    ASSERT_EQ(4u, ts.steps());
    ASSERT_EQ(8., ts.previous());
    ASSERT_EQ(12., ts.current());
    ASSERT_EQ(4., ts.dt());
    ASSERT_TRUE(alg.accepted());

    // dt*=0.5
    alg.setIterationNumber(7);
    ASSERT_TRUE(alg.next(solution_error));  // t=14, dt=2
    ts = alg.getTimeStep();
    ASSERT_EQ(5u, ts.steps());
    ASSERT_EQ(12., ts.previous());
    ASSERT_EQ(14., ts.current());
    ASSERT_EQ(2., ts.dt());
    ASSERT_TRUE(alg.accepted());

    // dt*=0.25 but dt_min = 1
    alg.setIterationNumber(8);              // exceed max
    ASSERT_TRUE(alg.next(solution_error));  // t=13, dt=1
    ts = alg.getTimeStep();
    ASSERT_EQ(5u, ts.steps());
    ASSERT_EQ(12., ts.previous());
    ASSERT_EQ(13, ts.current());
    ASSERT_EQ(1., ts.dt());
    ASSERT_FALSE(alg.accepted());

    // restart, dt*=1
    alg.setIterationNumber(4);
    ASSERT_TRUE(alg.next(solution_error));  // t=14, dt=1
    ts = alg.getTimeStep();
    ASSERT_EQ(6u, ts.steps());
    ASSERT_EQ(13., ts.previous());
    ASSERT_EQ(14, ts.current());
    ASSERT_EQ(1., ts.dt());
    ASSERT_TRUE(alg.accepted());
}

TEST(NumLib, TimeSteppingIterationNumberBased2)
{
    std::vector<std::size_t> iter_times_vector = {0, 3, 5, 7};
    std::vector<double> multiplier_vector = {2.0, 1.0, 0.5, 0.25};
    NumLib::IterationNumberBasedAdaptiveTimeStepping alg(1, 31, 1, 10, 1, iter_times_vector, multiplier_vector);

    std::vector<std::size_t> nr_iterations = {2, 2, 2, 4, 6, 8, 4, 4, 2, 2};
    const std::vector<double> expected_vec_t = {1, 2, 4, 8, 16, 24, 26, 28, 30, 31};

    struct IterationNumberUpdate
    {
        IterationNumberUpdate(std::vector<std::size_t> vec,
                              std::size_t& counter)
            : _nr_iterations(std::move(vec)), i(counter)
        {
        }

        std::vector<std::size_t> _nr_iterations;
        std::size_t& i;

        void operator()(NumLib::IterationNumberBasedAdaptiveTimeStepping &obj)
        {
            std::size_t n = (i<_nr_iterations.size()) ? _nr_iterations[i++] : 0;
            //INFO("-> NR-iterations=%d", n);
            obj.setIterationNumber(n);
        }
    };

    std::size_t counter = 0;
    IterationNumberUpdate update(nr_iterations, counter);

    std::vector<double> vec_t = timeStepping(alg, &update);
    //std::cout << vec_t;

    ASSERT_EQ(expected_vec_t.size(), vec_t.size());
    ASSERT_EQ(1u, alg.getNumberOfRepeatedSteps());
    ASSERT_ARRAY_NEAR(expected_vec_t, vec_t, expected_vec_t.size(), std::numeric_limits<double>::epsilon());
}
