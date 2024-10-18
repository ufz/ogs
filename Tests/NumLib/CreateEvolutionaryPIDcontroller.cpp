/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

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
