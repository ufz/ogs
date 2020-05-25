/**
 * \author Norihiro Watanabe
 * \date   2012-08-03
 *
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
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
        : previous_(current_time),
          current_(current_time),
          dt_(current_ - previous_),
          steps_(0)
    {
    }

    /**
     * Initialize a time step
     * @param previous_time    previous time
     * @param current_time     current time
     * @param n                the number of time steps
     */
    TimeStep(double previous_time, double current_time, std::size_t n)
        : previous_(previous_time),
          current_(current_time),
          dt_(current_ - previous_),
          steps_(n)
    {
    }

    /// copy a time step
    TimeStep(const TimeStep& src)
        : previous_(src.previous_),
          current_(src.current_),
          dt_(current_ - previous_),
          steps_(src.steps_)
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
        previous_ = current_;
        current_ += dt;
        dt_ = dt;
        steps_++;
        return *this;
    }

    /// compare current time
    bool operator<(const TimeStep& t) const { return (current_ < t.current_); }
    /// compare current time
    bool operator<(const double& t) const { return (current_ < t); }
    /// compare current time
    bool operator<=(const TimeStep& t) const
    {
        return (current_ <= t.current_);
    }

    /// compare current time
    bool operator<=(const double& t) const { return (current_ <= t); }
    /// compare current time
    bool operator==(const TimeStep& t) const
    {
        return (current_ == t.current_);
    }

    /// compare current time
    bool operator==(const double& t) const { return (current_ == t); }
    /// return previous time step
    double previous() const { return previous_; }
    /// return current time step
    double current() const { return current_; }
    /// time step size from previous_
    double dt() const { return dt_; }
    /// the number of time steps_
    std::size_t steps() const { return steps_; }

private:
    /// previous time step
    double previous_;
    /// current time step
    double current_;
    /// time step size
    double dt_;
    /// the number of time steps
    std::size_t steps_;
};

}  // namespace NumLib
