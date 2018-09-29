/**
 * \author Norihiro Watanabe
 * \date   2012-08-03
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <gtest/gtest.h>

#include <vector>

#include <logog/include/logog.hpp>

#include "NumLib/TimeStepping/TimeStep.h"
#include "NumLib/TimeStepping/Algorithms/FixedTimeStepping.h"

#include "Tests/TestTools.h"
#include "TimeSteppingTestingTools.h"

TEST(NumLib, TimeSteppingFixed)
{
    // homogeneous dt
    {
        NumLib::FixedTimeStepping fixed(1, 31, 10);
        const std::vector<double> expected_vec_t = {1, 11, 21, 31};

        std::vector<double> vec_t = timeStepping(fixed);

        ASSERT_EQ(expected_vec_t.size(), vec_t.size());
        ASSERT_ARRAY_NEAR(expected_vec_t, vec_t, expected_vec_t.size(), std::numeric_limits<double>::epsilon());
    }

    // dt vector (t_end == t0 + sum(dt))
    {
        const std::vector<double> fixed_dt = {10, 10, 10};
        NumLib::FixedTimeStepping fixed(1, 31, 10);
        const std::vector<double> expected_vec_t = {1, 11, 21, 31};

        std::vector<double> vec_t = timeStepping(fixed);

        ASSERT_EQ(expected_vec_t.size(), vec_t.size());
        ASSERT_ARRAY_NEAR(expected_vec_t, vec_t, expected_vec_t.size(), std::numeric_limits<double>::epsilon());
    }

    // dt vector (t_end < t0 + sum(dt))
    {
        const std::vector<double> fixed_dt = {5, 10, 20};
        NumLib::FixedTimeStepping fixed(1, 31, fixed_dt);
        const std::vector<double> expected_vec_t = {1, 6, 16, 31};

        std::vector<double> vec_t = timeStepping(fixed);

        ASSERT_EQ(expected_vec_t.size(), vec_t.size());
        ASSERT_ARRAY_NEAR(expected_vec_t, vec_t, expected_vec_t.size(), std::numeric_limits<double>::epsilon());
    }

    // dt vector (t_end > t0 + sum(dt))
    {
        const std::vector<double> fixed_dt = {5, 10, 10};
        NumLib::FixedTimeStepping fixed(1, 31, fixed_dt);
        const std::vector<double> expected_vec_t = {1, 6, 16, 26};

        std::vector<double> vec_t = timeStepping(fixed);

        ASSERT_EQ(expected_vec_t.size(), vec_t.size());
        ASSERT_ARRAY_NEAR(expected_vec_t, vec_t, expected_vec_t.size(), std::numeric_limits<double>::epsilon());
    }
}
