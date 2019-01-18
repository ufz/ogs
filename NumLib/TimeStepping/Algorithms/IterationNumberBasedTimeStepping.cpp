/**
 * \file
 * \author Haibing Shao and Norihiro Watanabe
 * \date   2013-08-07
 *
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
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

bool IterationNumberBasedTimeStepping::next(double const /*solution_error*/,
                                            int const number_iterations)
{
    _iter_times = number_iterations;

    // check current time step
    if (std::abs(_ts_current.current() - end()) <
        std::numeric_limits<double>::epsilon())
    {
        return false;
    }

    // confirm current time and move to the next if accepted
    if (accepted())
    {
        _ts_prev = _ts_current;
        _dt_vector.push_back(_ts_current.dt());
    }
    else
    {
        ++_n_rejected_steps;
    }

    // prepare the next time step info
    _ts_current = _ts_prev;
    _ts_current += getNextTimeStepSize();

    return true;
}

double IterationNumberBasedTimeStepping::findMultiplier(
    int const number_iterations) const
{
    double multiplier = 1.0;
    // get the first multiplier by default
    if (!_multiplier_vector.empty())
    {
        multiplier = _multiplier_vector[0];
    }
    // finding the right multiplier
    for (std::size_t i = 0; i < _iter_times_vector.size(); i++)
    {
        if (_iter_times >= _iter_times_vector[i])
        {
            multiplier = _multiplier_vector[i];
        }
    }
    return multiplier;
}

double IterationNumberBasedTimeStepping::getNextTimeStepSize() const
{
    double dt = 0.0;

    // In first time step and first non-linear iteration take the initial dt.
    if (_ts_prev.steps() == 0 && _iter_times == 0)
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

    dt = std::min(std::max(dt, _min_dt), _max_dt);

    double const t_next = dt + _ts_prev.current();
    if (t_next > end())
    {
        dt = end() - _ts_prev.current();
    }

    return dt;
}

bool IterationNumberBasedTimeStepping::accepted() const
{
    return _iter_times <= _max_iter;
}
}  // namespace NumLib
