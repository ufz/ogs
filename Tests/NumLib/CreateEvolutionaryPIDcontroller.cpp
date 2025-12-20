// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "NumLib/TimeStepping/Algorithms/CreateEvolutionaryPIDcontroller.h"

#include <gtest/gtest.h>

#include <numeric>

#include "NumLib/TimeStepping/Algorithms/TimeStepAlgorithm.h"

TEST(NumLibCreateEvolutionaryPIDcontrollerTimeStepping,
     InitialTimeGreaterThanEndTime)
{
    double const t_initial = 1;
    double const t_end = 0;

    NumLib::EvolutionaryPIDcontrollerParameters parameters{
        t_initial, t_end, 0, 0, 0, 0, 0, 0};
    EXPECT_ANY_THROW(
        auto fixed_time_step_algorithm =
            NumLib::createEvolutionaryPIDcontroller(std::move(parameters), {}));
}
