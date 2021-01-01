/**
 * \author Norihiro Watanabe
 * \date   2012-08-03
 *
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
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

std::tuple<bool, double> FixedTimeStepping::next(
    double const /*solution_error*/, int const /*number_iterations*/)
{
    // check if last time step
    if (_ts_current.steps() == _dt_vector.size() ||
        std::abs(_ts_current.current() - _t_end) <
            std::numeric_limits<double>::epsilon())
    {
        return std::make_tuple(false, 0.0);
    }

    double dt = _dt_vector[_ts_current.steps()];
    if (_ts_current.current() + dt > _t_end)
    {  // upper bound by t_end
        dt = _t_end - _ts_current.current();
    }

    return std::make_tuple(true, dt);
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
