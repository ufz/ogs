/**
 * \author Haibing Shao and Norihiro Watanabe
 * \date   2013-08-07
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "IterationNumberBasedAdaptiveTimeStepping.h"

#include <limits>
#include <algorithm>
#include <cassert>
#include <cmath>

namespace NumLib
{

IterationNumberBasedAdaptiveTimeStepping::IterationNumberBasedAdaptiveTimeStepping(double t0, double tn,
                                double min_ts, double max_ts, double initial_ts,
                                const std::vector<std::size_t> &iter_times_vector,
                                const std::vector<double> &multiplier_vector)
: _t_initial(t0), _t_end(tn), _iter_times_vector(iter_times_vector), _multiplier_vector(multiplier_vector),
  _min_ts(min_ts), _max_ts(max_ts), _initial_ts(initial_ts), _max_iter(_iter_times_vector.empty() ? 0 : _iter_times_vector.back()),
  _iter_times(0), _ts_pre(t0), _ts_current(t0), _n_rejected_steps(0)
{
    assert(iter_times_vector.size() == multiplier_vector.size());
}

bool IterationNumberBasedAdaptiveTimeStepping::next()
{
    // check current time step
    if (std::abs(_ts_current.current()-_t_end) < std::numeric_limits<double>::epsilon())
        return false;

    // confirm current time and move to the next if accepted
    if (accepted()) {
        _ts_pre = _ts_current;
        _dt_vector.push_back(_ts_current.dt());
    } else {
        ++_n_rejected_steps;
    }

    // prepare the next time step info
    _ts_current = _ts_pre;
    _ts_current += getNextTimeStepSize();

    return true;
}

const TimeStep IterationNumberBasedAdaptiveTimeStepping::getTimeStep() const
{
    return _ts_current;
}

double IterationNumberBasedAdaptiveTimeStepping::getNextTimeStepSize() const
{
    double dt = 0.0;

    // if this is the first time step
    // then we use initial guess provided by a user
    if ( _ts_pre.steps() == 0 )
    {
        dt = _initial_ts;
    }
    else  // not the first time step
    {
        double tmp_multiplier = 1.0;
        // get the first multiplier by default
        if (!_multiplier_vector.empty())
            tmp_multiplier = _multiplier_vector[0];
        // finding the right multiplier
        for (std::size_t i=0; i<_iter_times_vector.size(); i++ )
            if ( this->_iter_times > _iter_times_vector[i] )
                tmp_multiplier = _multiplier_vector[i];
        // multiply the the multiplier
        dt = _ts_pre.dt() * tmp_multiplier;
    }

    // check whether out of the boundary
    if ( dt < _min_ts )
        dt = _min_ts;
    else if ( dt > _max_ts )
        dt = _max_ts;

    double t_next = dt + _ts_pre.current();
    if (t_next > end())
        dt = end() - _ts_pre.current();

    return dt;
}

bool IterationNumberBasedAdaptiveTimeStepping::accepted() const
{
    return ( this->_iter_times <= this->_max_iter );
}

} // NumLib
