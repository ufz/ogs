// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "NumLib/TimeStepping/Algorithms/CreateIterationNumberBasedTimeStepping.h"

#include <gtest/gtest.h>

#include <numeric>

#include "NumLib/TimeStepping/Algorithms/MultiplyerInterpolationType.h"
#include "NumLib/TimeStepping/Algorithms/TimeStepAlgorithm.h"

// No time steps can be generated for an empty time interval.
struct NumLibCreateIterationNumberBasedTimeSteppingEmptyTimeInterval
    : public ::testing::TestWithParam<std::pair<double, double>>
{
};

TEST_P(NumLibCreateIterationNumberBasedTimeSteppingEmptyTimeInterval,
       CreationThrows)
{
    auto const [t_initial, t_end] = GetParam();

    NumLib::IterationNumberBasedTimeSteppingParameters parameters{
        t_initial, t_end,
        0,         0,
        0,         NumLib::MultiplyerInterpolationType::PiecewiseConstant,
        {},        {}};
    EXPECT_ANY_THROW(auto fixed_time_step_algorithm =
                         NumLib::createIterationNumberBasedTimeStepping(
                             std::move(parameters), {}));
}

INSTANTIATE_TEST_SUITE_P(
    NumLibCreateIterationNumberBasedTimeStepping,
    NumLibCreateIterationNumberBasedTimeSteppingEmptyTimeInterval,
    // Initial time greater than the end time, initial time equal to it.
    ::testing::Values(std::pair{1.0, 0.0}, std::pair{1.0, 1.0}));
