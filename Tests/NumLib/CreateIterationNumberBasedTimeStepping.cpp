// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "NumLib/TimeStepping/Algorithms/CreateIterationNumberBasedTimeStepping.h"

#include <gtest/gtest.h>

#include <numeric>

#include "NumLib/TimeStepping/Algorithms/MultiplyerInterpolationType.h"
#include "NumLib/TimeStepping/Algorithms/TimeStepAlgorithm.h"

TEST(NumLibCreateIterationNumberBasedTimeStepping,
     InitialTimeGreaterThanEndTime)
{
    double const t_initial = 1;
    double const t_end = 0;

    NumLib::IterationNumberBasedTimeSteppingParameters parameters{
        t_initial, t_end,
        0,         0,
        0,         NumLib::MultiplyerInterpolationType::PiecewiseConstant,
        {},        {}};
    EXPECT_ANY_THROW(auto fixed_time_step_algorithm =
                         NumLib::createIterationNumberBasedTimeStepping(
                             std::move(parameters), {}));
}
