/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "TimeStepAlgorithm.h"

#include <algorithm>
#include <limits>

#include "NumLib/TimeStepping/Time.h"

namespace NumLib
{
double possiblyClampDtToNextFixedTime(
    Time const& t, double const dt,
    std::vector<double> const& fixed_output_times)
{
    auto const specific_time = std::upper_bound(
        std::cbegin(fixed_output_times), std::cend(fixed_output_times), t());

    if (specific_time == std::cend(fixed_output_times))
    {
        return dt;
    }

    if ((t < Time(*specific_time)) && t + dt > Time(*specific_time))
    {
        double const t_to_specific_time = *specific_time - t();
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
