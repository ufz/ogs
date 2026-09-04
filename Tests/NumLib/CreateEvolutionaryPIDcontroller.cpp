// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "NumLib/TimeStepping/Algorithms/CreateEvolutionaryPIDcontroller.h"

#include <gtest/gtest.h>

#include <numeric>

#include "NumLib/TimeStepping/Algorithms/TimeStepAlgorithm.h"

// No time steps can be generated for an empty time interval.
struct NumLibCreateEvolutionaryPIDcontrollerEmptyTimeInterval
    : public ::testing::TestWithParam<std::pair<double, double>>
{
};

TEST_P(NumLibCreateEvolutionaryPIDcontrollerEmptyTimeInterval, CreationThrows)
{
    auto const [t_initial, t_end] = GetParam();

    NumLib::EvolutionaryPIDcontrollerParameters parameters{
        t_initial, t_end, 0, 0, 0, 0, 0, 0};
    EXPECT_ANY_THROW(
        auto fixed_time_step_algorithm =
            NumLib::createEvolutionaryPIDcontroller(std::move(parameters), {}));
}

INSTANTIATE_TEST_SUITE_P(
    NumLibCreateEvolutionaryPIDcontrollerTimeStepping,
    NumLibCreateEvolutionaryPIDcontrollerEmptyTimeInterval,
    // Initial time greater than the end time, initial time equal to it.
    ::testing::Values(std::pair{1.0, 0.0}, std::pair{1.0, 1.0}));
