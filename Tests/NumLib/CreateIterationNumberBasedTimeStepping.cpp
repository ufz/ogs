/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "NumLib/TimeStepping/Algorithms/CreateIterationNumberBasedTimeStepping.h"

#include <gtest/gtest.h>

#include <numeric>

#include "NumLib/TimeStepping/Algorithms/TimeStepAlgorithm.h"

TEST(NumLibCreateIterationNumberBasedTimeStepping,
     InitialTimeGreaterThanEndTime)
{
    double const t_initial = 1;
    double const t_end = 0;

    NumLib::IterationNumberBasedTimeSteppingParameters parameters{
        t_initial, t_end, 0, 0, 0, {}, {}};
    EXPECT_ANY_THROW(auto fixed_time_step_algorithm =
                         NumLib::createIterationNumberBasedTimeStepping(
                             std::move(parameters), {}));
}
