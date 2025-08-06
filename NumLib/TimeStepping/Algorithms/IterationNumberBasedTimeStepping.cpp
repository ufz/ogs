/**
 * \file
 * \author Haibing Shao and Norihiro Watanabe
 * \date   2013-08-07
 *
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "IterationNumberBasedTimeStepping.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <utility>

#include "BaseLib/Algorithm.h"
#include "IterationNumberBasedTimeStepping.h"
#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"

namespace NumLib
{
IterationNumberBasedTimeStepping::IterationNumberBasedTimeStepping(
    double const t_initial, double const t_end, double const min_dt,
    double const max_dt, double const initial_dt,
    MultiplyerInterpolationType const multiplier_interpolation_type,
    std::vector<int>&& iter_times_vector,
    std::vector<double>&& multiplier_vector,
    std::vector<double> const& fixed_times_for_output)
    : TimeStepAlgorithm(t_initial, t_end),
      _iter_times_vector(std::move(iter_times_vector)),
      _multiplier_vector(std::move(multiplier_vector)),
      _min_dt(min_dt),
      _max_dt(max_dt),
      _initial_dt(initial_dt),
      _multiplier_interpolation_type(multiplier_interpolation_type),
      _max_iter(_iter_times_vector.empty() ? 0 : _iter_times_vector.back()),
      _fixed_times_for_output(fixed_times_for_output)
{
    if (_iter_times_vector.empty())
    {
        OGS_FATAL("Vector of iteration numbers must not be empty.");
    }
    if (_iter_times_vector.size() != _multiplier_vector.size())
    {
        OGS_FATAL(
            "Vector of iteration numbers must be of the same size as the "
            "vector of multipliers.");
    }
    if (!std::is_sorted(std::begin(_iter_times_vector),
                        std::end(_iter_times_vector)))
    {
        OGS_FATAL("Vector of iteration numbers must be sorted.");
    }
}

double IterationNumberBasedTimeStepping::next(double const /*solution_error*/,
                                              int const number_iterations,
                                              NumLib::TimeStep& ts_previous,
                                              NumLib::TimeStep& ts_current)
{
    _iter_times = number_iterations;

    if (ts_previous.isAccepted())
    {
        ts_previous = ts_current;
    }

    // confirm current time and move to the next if accepted
    if (ts_current.isAccepted())
    {
        ts_previous.setAccepted(true);
        return getNextTimeStepSize(ts_previous, ts_current);
    }

    double dt = getNextTimeStepSize(ts_previous, ts_current);
    // In case it is the first time be rejected, re-computed dt again with
    // current dt
    if (std::abs(dt - ts_current.dt()) < std::numeric_limits<double>::epsilon())
    {
        // time step was rejected, keep dt for the next dt computation.
        ts_previous =  // essentially equal to _ts_prev.dt = _ts_current.dt.
            TimeStep{ts_previous.previous(), ts_previous.previous() + dt,
                     ts_previous.timeStepNumber()};
        dt = getNextTimeStepSize(ts_previous, ts_current);
    }

    // time step was rejected, keep dt for the next dt computation.
    ts_previous =  // essentially equal to ts_previous.dt = _ts_current.dt.
        TimeStep{ts_previous.previous(), ts_previous.previous() + dt,
                 ts_previous.timeStepNumber()};
    ts_current = TimeStep{ts_current.previous(), ts_current.previous() + dt,
                          ts_current.timeStepNumber()};

    return dt;
}

double findMultiplier(
    int const number_iterations, bool const current_time_step_is_accepted,
    std::vector<int> const& nonlinear_iteration_numbers,
    std::vector<double> const& multipliers,
    MultiplyerInterpolationType const multiplier_interpolation_type)
{
    double multiplier = multipliers.front();
    switch (multiplier_interpolation_type)
    {
        case MultiplyerInterpolationType::PiecewiseLinear:
        {
            auto const& pwli = MathLib::PiecewiseLinearInterpolation(
                nonlinear_iteration_numbers, multipliers, false);
            multiplier = pwli.getValue(number_iterations);
            DBUG("Using piecewise linear iteration-based time stepping.");
            break;
        }
        case MultiplyerInterpolationType::PiecewiseConstant:
            DBUG("Using piecewise constant iteration-based time stepping.");
            for (std::size_t i = 0; i < nonlinear_iteration_numbers.size(); i++)
            {
                if (number_iterations >= nonlinear_iteration_numbers[i])
                {
                    multiplier = multipliers[i];
                }
            }
            break;
    }

    if (!current_time_step_is_accepted && (multiplier >= 1.0))
    {
        return *std::min_element(multipliers.begin(), multipliers.end());
    }

    return multiplier;
}

double IterationNumberBasedTimeStepping::getNextTimeStepSize(
    NumLib::TimeStep const& ts_previous,
    NumLib::TimeStep const& ts_current) const
{
    double dt = 0.0;

    // In first time step and first non-linear iteration take the initial dt.
    if (ts_previous.timeStepNumber() == 0 && _iter_times == 0)
    {
        dt = _initial_dt;
    }
    else
    {
        // Attention: for the first time step and second iteration the
        // ts_prev.dt is 0 and 0*multiplier is the next dt, which will be
        // clamped to the minimum dt.
        dt = ts_previous.dt() *
             findMultiplier(_iter_times, ts_current.isAccepted(),
                            _iter_times_vector, _multiplier_vector,
                            _multiplier_interpolation_type);
    }

    if (_fixed_times_for_output.empty())
    {
        return std::clamp(dt, _min_dt, _max_dt);
    }

    // restrict dt to _max_dt before taking fixed times for output into account
    dt = std::min(dt, _max_dt);

    // find first fixed timestep for output larger than the current time, i.e.,
    // current time < fixed output time
    auto fixed_output_time_it = std::find_if(
        std::begin(_fixed_times_for_output), std::end(_fixed_times_for_output),
        [&ts_current](auto const fixed_output_time)
        { return ts_current.current()() < fixed_output_time; });

    if (fixed_output_time_it != _fixed_times_for_output.end())
    {
        // check if the fixed output time is in the interval
        // (current time, current time + dt)
        if (*fixed_output_time_it < ts_current.current()() + dt)
        {
            // check if the potential adjusted time step is larger than zero
            if (std::abs(*fixed_output_time_it - ts_current.current()()) >
                std::numeric_limits<double>::epsilon() * ts_current.current()())
            {
                return *fixed_output_time_it - ts_current.current()();
            }
        }
    }
    return std::clamp(dt, _min_dt, _max_dt);
}

bool IterationNumberBasedTimeStepping::canReduceTimestepSize(
    NumLib::TimeStep const& timestep_previous,
    NumLib::TimeStep const& timestep_current) const
{
    return NumLib::canReduceTimestepSize(timestep_previous, timestep_current,
                                         _min_dt);
}

}  // namespace NumLib
