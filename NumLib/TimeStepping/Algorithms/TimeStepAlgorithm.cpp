/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "TimeStepAlgorithm.h"

#include <algorithm>

namespace NumLib
{
double possiblyClampDtToNextFixedTime(
    double const t, double const dt,
    std::vector<double> const& fixed_output_times)
{
    auto const specific_time = std::upper_bound(
        std::cbegin(fixed_output_times), std::cend(fixed_output_times), t);

    if (specific_time == std::cend(fixed_output_times))
    {
        return dt;
    }

    if ((*specific_time > t) && (t + dt - *specific_time > 0.0))
    {
        return *specific_time - t;
    }

    return dt;
}
}  // namespace NumLib
