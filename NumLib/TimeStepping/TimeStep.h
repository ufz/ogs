/**
 * \author Norihiro Watanabe
 * \date   2012-08-03
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
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
class TimeStep
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
          _steps(0)
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
          _steps(n)
    {
    }

    /// copy a time step
    TimeStep(const TimeStep& src)
        : _previous(src._previous),
          _current(src._current),
          _dt(_current - _previous),
          _steps(src._steps)
    {
    }

    /// copy a time step
    TimeStep& operator=(const TimeStep& src) = default;

    /// return a time step incremented by the given time step size
    TimeStep operator+(const double dt) const
    {
        TimeStep t(*this);
        t += dt;
        return t;
    }

    /// increment time step
    TimeStep& operator+=(const double dt)
    {
        _previous = _current;
        _current += dt;
        _dt = dt;
        _steps++;
        return *this;
    }

    /// compare current time
    bool operator<(const TimeStep& t) const { return (_current < t._current); }
    /// compare current time
    bool operator<(const double& t) const { return (_current < t); }
    /// compare current time
    bool operator<=(const TimeStep& t) const
    {
        return (_current <= t._current);
    }

    /// compare current time
    bool operator<=(const double& t) const { return (_current <= t); }
    /// compare current time
    bool operator==(const TimeStep& t) const
    {
        return (_current == t._current);
    }

    /// compare current time
    bool operator==(const double& t) const { return (_current == t); }
    /// return previous time step
    double previous() const { return _previous; }
    /// return current time step
    double current() const { return _current; }
    /// time step size from _previous
    double dt() const { return _dt; }
    /// the number of time _steps
    std::size_t steps() const { return _steps; }

private:
    /// previous time step
    double _previous;
    /// current time step
    double _current;
    /// time step size
    double _dt;
    /// the number of time steps
    std::size_t _steps;
};

}  // NumLib
