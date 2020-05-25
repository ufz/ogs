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

#include "FixedTimeStepping.h"

#include <algorithm>
#include <numeric>
#include <limits>
#include <cassert>

namespace NumLib
{
FixedTimeStepping::FixedTimeStepping(double t0,
                                     double tn,
                                     const std::vector<double>& vec_all_dt)
    : TimeStepAlgorithm(t0, computeEnd(t0, tn, vec_all_dt), vec_all_dt)
{
}

FixedTimeStepping::FixedTimeStepping(double t0, double tn, double dt)
    : TimeStepAlgorithm(t0, tn, dt)
{
}

bool FixedTimeStepping::next(double const /*solution_error*/,
                             int const /*number_iterations*/)
{
    // check if last time step
    if (ts_current_.steps() == dt_vector_.size() ||
        std::abs(ts_current_.current() - t_end_) <
            std::numeric_limits<double>::epsilon())
    {
        return false;
    }

    // confirm current time and move to the next if accepted
    if (accepted())
    {
        ts_prev_ = ts_current_;
    }

    // prepare the next time step info
    ts_current_ = ts_prev_;
    double dt = dt_vector_[ts_prev_.steps()];
    if (ts_prev_.current() + dt > t_end_)
    {  // upper bound by t_end
        dt = t_end_ - ts_prev_.current();
    }
    ts_current_ += dt;

    return true;
}

double FixedTimeStepping::computeEnd(double t_initial,
                                     double t_end,
                                     const std::vector<double>& dt_vector)
{
    double t_sum =
        t_initial + std::accumulate(dt_vector.begin(), dt_vector.end(), 0.);
    return std::min(t_end, t_sum);
}

}  // namespace NumLib
