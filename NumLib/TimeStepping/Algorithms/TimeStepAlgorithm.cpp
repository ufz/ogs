/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "TimeStepAlgorithm.h"

#include <range/v3/algorithm/upper_bound.hpp>

#include "NumLib/TimeStepping/Time.h"

namespace NumLib
{
double possiblyClampDtToNextFixedTime(
    Time const& t, double const dt,
    std::vector<double> const& fixed_output_times)
{
    auto const specific_time = ranges::upper_bound(
        fixed_output_times, t, ranges::less{}, [](auto t) { return Time(t); });

    if (specific_time == ranges::cend(fixed_output_times))
    {
        return dt;
    }

    Time const fixed_output_time(*specific_time);
    if ((t < fixed_output_time) && (t + dt) > fixed_output_time)
    {
        double const t_to_specific_time = fixed_output_time() - t();
        return t_to_specific_time;
    }

    return dt;
}

bool canReduceTimestepSize(TimeStep const& timestep_previous,
                           TimeStep const& timestep_current,
                           double const min_dt)
{
    return !(timestep_current.dt() == min_dt &&
             timestep_previous.dt() == min_dt);
}
}  // namespace NumLib
