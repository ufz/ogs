/**
 * \author Norihiro Watanabe
 * \date   2012-08-03
 *
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <cstddef>

namespace NumLib
{
/**
 * \brief Time step object
 *
 * TimeStep class contains previous time(\f$t_{n}\f$), current time
 * (\f$t_{n+1}\f$),
 * time step length (\f$\Delta t_{n+1}\f$), and the number of time steps
 * (\f$n+1\f$).
 */
class TimeStep final
{
public:
    /**
     * Initialize a time step
     * @param current_time     current time
     */
    explicit TimeStep(double current_time)
        : _previous(current_time),
          _current(current_time),
          _dt(_current - _previous),
          _time_step_number(0)
    {
    }

    /**
     * Initialize a time step
     * @param previous_time    previous time
     * @param current_time     current time
     * @param n                the number of time steps
     */
    TimeStep(double previous_time, double current_time, std::size_t n)
        : _previous(previous_time),
          _current(current_time),
          _dt(_current - _previous),
          _time_step_number(n)
    {
    }

    /// copy a time step
    TimeStep(const TimeStep& src) = default;

    /// copy a time step
    TimeStep& operator=(const TimeStep& src) = default;

    /// increment time step
    TimeStep& operator+=(const double dt)
    {
        _previous = _current;
        _current += dt;
        _dt = dt;
        _time_step_number++;
        return *this;
    }

    /// compare current time
    bool operator<(TimeStep const& ts) const
    {
        return (_current < ts._current);
    }
    /// compare current time
    bool operator<=(TimeStep const& ts) const
    {
        return (_current <= ts._current);
    }
    /// compare current time
    bool operator==(TimeStep const& ts) const
    {
        return (_current == ts._current);
    }
    /// return previous time step
    double previous() const { return _previous; }
    /// return current time step
    double current() const { return _current; }
    /// time step size from _previous
    double dt() const { return _dt; }
    /// the time step number
    std::size_t timeStepNumber() const { return _time_step_number; }

private:
    /// previous time step
    double _previous;
    /// current time step
    double _current;
    /// time step size
    double _dt;
    /// the number of time steps
    std::size_t _time_step_number;
};

}  // namespace NumLib
