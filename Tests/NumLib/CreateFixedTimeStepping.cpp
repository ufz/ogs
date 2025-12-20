// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "NumLib/TimeStepping/Algorithms/CreateFixedTimeStepping.h"

#include <gtest/gtest.h>

#include <numeric>

TEST(NumLibCreateFixedTimeStepping, InitialTimeGreaterThanEndTime)
{
    double const t_initial = 1;
    double const t_end = 0;
    std::vector<NumLib::RepeatDtPair> const repeat_dt_pair;

    NumLib::FixedTimeSteppingParameters const parameters{t_initial, t_end,
                                                         repeat_dt_pair};
    EXPECT_ANY_THROW(auto fixed_time_step_algorithm =
                         NumLib::createFixedTimeStepping(parameters, {}));
}
