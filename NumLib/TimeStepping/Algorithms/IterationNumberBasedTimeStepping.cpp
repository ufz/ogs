/**
 * \file
 * \author Haibing Shao and Norihiro Watanabe
 * \date   2013-08-07
 *
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
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

namespace NumLib
{
IterationNumberBasedTimeStepping::IterationNumberBasedTimeStepping(
    double const t_initial, double const t_end, double const min_dt,
    double const max_dt, double const initial_dt,
    std::vector<int>&& iter_times_vector,
    std::vector<double>&& multiplier_vector)
    : TimeStepAlgorithm(t_initial, t_end),
      _iter_times_vector(std::move(iter_times_vector)),
      _multiplier_vector(std::move(multiplier_vector)),
      _min_dt(min_dt),
      _max_dt(max_dt),
      _initial_dt(initial_dt),
      _max_iter(_iter_times_vector.empty() ? 0 : _iter_times_vector.back())
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

std::tuple<bool, double> IterationNumberBasedTimeStepping::next(
    double const /*solution_error*/, int const number_iterations)
{
    _iter_times = number_iterations;

    if (_previous_time_step_accepted)
    {
        _ts_prev = _ts_current;
    }

    // confirm current time and move to the next if accepted
    if (accepted())
    {
        _previous_time_step_accepted = true;
        return std::make_tuple(true, getNextTimeStepSize());
    }
    else
    {
        double dt = getNextTimeStepSize();
        // In case it is the first time be rejected, re-computed dt again with
        // current dt
        if (std::abs(dt - _ts_current.dt()) <
            std::numeric_limits<double>::epsilon())
        {
            // time step was rejected, keep dt for the next dt computation.
            _ts_prev =  // essentially equal to _ts_prev.dt = _ts_current.dt.
                TimeStep{_ts_prev.previous(), _ts_prev.previous() + dt,
                         _ts_prev.timeStepNumber()};
            dt = getNextTimeStepSize();
        }

        // time step was rejected, keep dt for the next dt computation.
        _ts_prev =  // essentially equal to _ts_prev.dt = _ts_current.dt.
            TimeStep{_ts_prev.previous(), _ts_prev.previous() + dt,
                     _ts_prev.timeStepNumber()};

        _previous_time_step_accepted = false;

        return std::make_tuple(false, dt);
    }
    return {};
}

double IterationNumberBasedTimeStepping::findMultiplier(
    int const number_iterations) const
{
    double multiplier = _multiplier_vector.front();
    for (std::size_t i = 0; i < _iter_times_vector.size(); i++)
    {
        if (number_iterations >= _iter_times_vector[i])
        {
            multiplier = _multiplier_vector[i];
        }
    }

    if (!_is_accepted && (multiplier >= 1.0))
    {
        return *std::min_element(_multiplier_vector.begin(),
                                 _multiplier_vector.end());
    }

    return multiplier;
}

double IterationNumberBasedTimeStepping::getNextTimeStepSize() const
{
    double dt = 0.0;

    // In first time step and first non-linear iteration take the initial dt.
    if (_ts_prev.timeStepNumber() == 0 && _iter_times == 0)
    {
        dt = _initial_dt;
    }
    else
    {
        // Attention: for the first time step and second iteration the
        // ts_prev.dt is 0 and 0*multiplier is the next dt, which will be
        // clamped to the minimum dt.
        dt = _ts_prev.dt() * findMultiplier(_iter_times);
    }

    return std::clamp(dt, _min_dt, _max_dt);
}

bool IterationNumberBasedTimeStepping::canReduceTimestepSize() const
{
    // If current and previous dt are both at minimum dt, then cannot reduce
    // further.
    return !(_ts_current.dt() == _min_dt && _ts_prev.dt() == _min_dt);
}

}  // namespace NumLib
