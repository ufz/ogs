/**
 * \file
 * \author Norihiro Watanabe
 * \date   2012-08-03
 *
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <tuple>
#include <vector>

#include "BaseLib/Logging.h"
#include "NumLib/TimeStepping/Algorithms/TimeStepAlgorithm.h"
#include "NumLib/TimeStepping/Time.h"
#include "NumLib/TimeStepping/TimeStep.h"

namespace
{

template <class T_TIME_STEPPING>
std::vector<double> timeStepping(T_TIME_STEPPING& algorithm,
                                 std::vector<int> const& number_iterations,
                                 std::vector<double> const& fixed_output_times)
{
    std::vector<double> vec_t;
    vec_t.push_back(algorithm.begin()());

    auto const end_time = algorithm.end();
    NumLib::TimeStep current_timestep(algorithm.begin());
    NumLib::TimeStep previous_timestep(algorithm.begin());

    double const solution_error = 0;
    for (auto const& i : number_iterations)
    {
        auto [step_accepted, timestepper_dt] = algorithm.next(
            solution_error, i, previous_timestep, current_timestep);
        if (!step_accepted)
        {
            break;
        }
        if (current_timestep.current() + timestepper_dt ==
            current_timestep.current())
        {
            break;
        }

        if (!fixed_output_times.empty())
        {
            timestepper_dt = NumLib::possiblyClampDtToNextFixedTime(
                current_timestep.current(), timestepper_dt, fixed_output_times);
        }

        timestepper_dt =
            (current_timestep.current() + timestepper_dt > end_time)
                ? end_time() - current_timestep.current()()
                : timestepper_dt;

        NumLib::updateTimeSteps(timestepper_dt, previous_timestep,
                                current_timestep);
        // INFO("t: n={:d},t={:g},dt={:g}", t.steps(), t.current(), t.dt());
        if (current_timestep.isAccepted())
        {
            vec_t.push_back(current_timestep.current()());
        }
        else
        {
            //INFO("*** rejected.");
        }
    }

    return vec_t;
}
} // namespace
