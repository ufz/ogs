/**
 * \file
 * \author Haibing Shao and Norihiro Watanabe
 * \date   2013-08-07
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
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
    std::vector<double>&& multiplier_vector,
    std::vector<double>&& fixed_output_times)
    : TimeStepAlgorithm(t_initial, t_end),
      _iter_times_vector(std::move(iter_times_vector)),
      _multiplier_vector(std::move(multiplier_vector)),
      _fixed_output_times(std::move(fixed_output_times)),
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

    // Remove possible duplicated elements. Result will be sorted.
    BaseLib::makeVectorUnique(_fixed_output_times);
}

bool IterationNumberBasedTimeStepping::next(double const /*solution_error*/,
                                            int const number_iterations)
{
    _iter_times = number_iterations;

    // confirm current time and move to the next if accepted
    if (accepted())
    {
        _ts_prev = _ts_current;
        _dt_vector.push_back(_ts_current.dt());
    }
    else
    {
        ++_n_rejected_steps;
        // time step was rejected, keep dt for the next dt computation.
        _ts_prev =  // essentially equal to _ts_prev.dt = _ts_current.dt.
            TimeStep{_ts_prev.previous(),
                     _ts_prev.previous() + _ts_current.dt(), _ts_prev.steps()};
    }

    // prepare the next time step info
    _ts_current = _ts_prev;
    _ts_current += possiblyClampDtToNextFixedTime(
        _ts_current.current(), getNextTimeStepSize(), _fixed_output_times);

    return true;
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

    if (!_accepted && (multiplier >= 1.0))
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

    dt = std::clamp(dt, _min_dt, _max_dt);

    double const t_next = dt + _ts_prev.current();
    if (t_next > end())
    {
        dt = end() - _ts_prev.current();
    }

    return dt;
}

void IterationNumberBasedTimeStepping::addFixedOutputTimes(
    std::vector<double> const& extra_fixed_output_times)
{
    _fixed_output_times.insert(_fixed_output_times.end(),
                               extra_fixed_output_times.begin(),
                               extra_fixed_output_times.end());

    // Remove possible duplicated elements. Result will be sorted.
    BaseLib::makeVectorUnique(_fixed_output_times);
}

bool IterationNumberBasedTimeStepping::accepted() const
{
    return _accepted;
}

bool IterationNumberBasedTimeStepping::canReduceTimestepSize() const
{
    // If current and previous dt are both at minimum dt, then cannot reduce
    // further.
    return !(_ts_current.dt() == _min_dt && _ts_prev.dt() == _min_dt);
}

}  // namespace NumLib
