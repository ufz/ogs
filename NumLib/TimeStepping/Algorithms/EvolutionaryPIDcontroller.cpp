/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  Created on March 31, 2017, 4:13 PM
 */

#include "EvolutionaryPIDcontroller.h"

#include <functional>
#include <limits>
#include <vector>

#include "BaseLib/Algorithm.h"
#include "BaseLib/Logging.h"

namespace NumLib
{
std::tuple<bool, double> EvolutionaryPIDcontroller::next(
    double const solution_error, int const /*number_iterations*/,
    NumLib::TimeStep& /*timestep_previous*/, NumLib::TimeStep& timestep_current)
{
    const bool is_previous_step_accepted = timestep_current.isAccepted();

    const double e_n = solution_error;
    const double zero_threshold = std::numeric_limits<double>::epsilon();
    // step rejected.
    if (e_n > _tol)
    {
        timestep_current.setAccepted(false);

        double h_new = (e_n > zero_threshold)
                           ? timestep_current.dt() * _tol / e_n
                           : 0.5 * timestep_current.dt();

        h_new =
            limitStepSize(h_new, is_previous_step_accepted, timestep_current);

        WARN(
            "This step is rejected due to the relative change from the"
            " solution of the previous\n"
            "\t time step to the current solution exceeds the given tolerance"
            " of {:g}.\n"
            "\t This time step will be repeated with a new time step size of"
            " {:g}\n"
            "\t or the simulation will be halted.",
            _tol, h_new);

        return std::make_tuple(timestep_current.isAccepted(), h_new);
    }

    // step accepted.
    timestep_current.setAccepted(true);

    if (timestep_current.timeStepNumber() == 0)
    {
        _e_n_minus1 = e_n;

        return std::make_tuple(timestep_current.isAccepted(), _h0);
    }
    else
    {
        const double h_n = timestep_current.dt();
        double h_new = h_n;

        if (e_n > zero_threshold)
        {
            if (_e_n_minus1 > zero_threshold)
            {
                if (_e_n_minus2 > zero_threshold)
                {
                    h_new = std::pow(_e_n_minus1 / e_n, _kP) *
                            std::pow(_tol / e_n, _kI) *
                            std::pow(
                                _e_n_minus1 * _e_n_minus1 / (e_n * _e_n_minus2),
                                _kD) *
                            h_n;
                }
                else
                {
                    h_new = std::pow(_e_n_minus1 / e_n, _kP) *
                            std::pow(_tol / e_n, _kI) * h_n;
                }
            }
            else
            {
                h_new = std::pow(_tol / e_n, _kI) * h_n;
            }
        }

        h_new =
            limitStepSize(h_new, is_previous_step_accepted, timestep_current);

        _e_n_minus2 = _e_n_minus1;
        _e_n_minus1 = e_n;

        return std::make_tuple(timestep_current.isAccepted(), h_new);
    }
}

double EvolutionaryPIDcontroller::limitStepSize(
    const double h_new, const bool previous_step_accepted,
    NumLib::TimeStep const& timestep_current) const
{
    const double h_n = timestep_current.dt();
    // Forced the computed time step size in the given range
    // (see the formulas in the documentation of the class)
    const double h_in_range = std::max(_h_min, std::min(h_new, _h_max));
    // Limit the step size change ratio.
    double limited_h =
        std::max(_rel_h_min * h_n, std::min(h_in_range, _rel_h_max * h_n));

    if (!previous_step_accepted)
    {
        // If the last time step was rejected and the new predicted time step
        // size is identical to that of the previous rejected step, the new
        // step size is then reduced by half.
        if (std::abs(limited_h - timestep_current.dt()) <
            std::numeric_limits<double>::min())
        {
            limited_h = std::max(_h_min, 0.5 * limited_h);
        }

        // If the last time step was rejected and the new predicted time step
        // size is larger than the step size of the rejected step, the new step
        // size takes the half of the size of the rejected step. This could
        // happen when a time step is rejected due to a diverged non-linear
        // solver. In such case, this algorithm may give a large time step size
        // by using the diverged solution.
        if (limited_h > timestep_current.dt())
        {
            limited_h = std::max(_h_min, 0.5 * timestep_current.dt());
        }
    }

    if (_fixed_times_for_output.empty())
    {
        return limited_h;
    }

    // find first fixed timestep for output larger than the current time, i.e.,
    // current time < fixed output time
    auto fixed_output_time_it = std::find_if(
        std::begin(_fixed_times_for_output), std::end(_fixed_times_for_output),
        [&timestep_current](auto const fixed_output_time)
        { return timestep_current.current() < fixed_output_time; });

    if (fixed_output_time_it != _fixed_times_for_output.end())
    {
        // check if the fixed output time is in the interval
        // (current time, current time + limited_h)
        if (*fixed_output_time_it < timestep_current.current() + limited_h)
        {
            // check if the potential adjusted time step is larger than zero
            if (std::abs(*fixed_output_time_it - timestep_current.current()) >
                std::numeric_limits<double>::epsilon() *
                    timestep_current.current())
            {
                return *fixed_output_time_it - timestep_current.current();
            }
        }
    }
    return limited_h;
}

bool EvolutionaryPIDcontroller::canReduceTimestepSize(
    NumLib::TimeStep const& timestep_previous,
    NumLib::TimeStep const& timestep_current) const
{
    return NumLib::canReduceTimestepSize(timestep_previous, timestep_current,
                                         _h_min);
}
}  // namespace NumLib
