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
                                 std::vector<double> const& fixed_output_times,
                                 std::vector<int> const& rejected_time_steps)
{
    std::vector<double> vec_t;
    vec_t.push_back(algorithm.begin()());

    auto const end_time = algorithm.end();
    NumLib::TimeStep current_timestep(algorithm.begin());
    NumLib::TimeStep previous_timestep(algorithm.begin());

    double const solution_error = 0;
    int time_step_counter = 0;
    std::size_t idx = 0;
    bool last_time_step_rejected = false;
    double timestepper_dt = 0.0;

    for (auto const& i : number_iterations)
    {
        if (idx < rejected_time_steps.size() &&
            time_step_counter == rejected_time_steps[idx])
        {
            current_timestep.setAccepted(false);
            idx++;
            timestepper_dt = algorithm.next(
                solution_error, i, previous_timestep, current_timestep);
            last_time_step_rejected = true;
        }
        else
        {
            timestepper_dt = algorithm.next(
                solution_error, i, previous_timestep, current_timestep);
            time_step_counter++;
            last_time_step_rejected = false;
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

        if (!last_time_step_rejected)
        {
            NumLib::updateTimeSteps(timestepper_dt, previous_timestep,
                                    current_timestep);
        }
        else
        {
            current_timestep = previous_timestep;
        }
        vec_t.push_back(current_timestep.current()());
    }

    return vec_t;
}
} // namespace
