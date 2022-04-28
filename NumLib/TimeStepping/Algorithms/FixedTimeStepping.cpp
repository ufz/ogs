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

#include "FixedTimeStepping.h"

#include <algorithm>
#include <cassert>
#include <limits>
#include <numeric>

namespace NumLib
{
/// determine true end time
double computeEnd(double t_initial,
                  double t_end,
                  const std::vector<double>& dt_vector)
{
    double t_sum =
        t_initial + std::accumulate(dt_vector.begin(), dt_vector.end(), 0.);
    return std::min(t_end, t_sum);
}

FixedTimeStepping::FixedTimeStepping(double t0,
                                     double tn,
                                     const std::vector<double>& vec_all_dt)
    : TimeStepAlgorithm(t0, computeEnd(t0, tn, vec_all_dt)),
      _dt_vector(vec_all_dt)
{
}

FixedTimeStepping::FixedTimeStepping(double t0, double t_end, double dt)
    : TimeStepAlgorithm(t0, t_end)
{
    auto const new_size =
        static_cast<std::size_t>(std::ceil((t_end - t0) / dt));
    try
    {
        _dt_vector = std::vector<double>(new_size, dt);
    }
    catch (std::length_error const& e)
    {
        OGS_FATAL(
            "Resize of the time steps vector failed for the requested new "
            "size {:d}. Probably there is not enough memory ({:g} GiB "
            "requested).\n"
            "Thrown exception: {:s}",
            new_size, new_size * sizeof(double) / 1024. / 1024. / 1024.,
            e.what());
    }
    catch (std::bad_alloc const& e)
    {
        OGS_FATAL(
            "Allocation of the time steps vector failed for the requested "
            "size {:d}. Probably there is not enough memory ({:d} GiB "
            "requested).\n"
            "Thrown exception: {:s}",
            new_size,
            new_size * sizeof(double) / 1024. / 1024. / 1024.,
            e.what());
    }
}

std::tuple<bool, double> FixedTimeStepping::next(
    double const /*solution_error*/, int const /*number_iterations*/,
    NumLib::TimeStep& /*ts_previous*/, NumLib::TimeStep& ts_current)
{
    // check if last time step
    if (ts_current.timeStepNumber() == _dt_vector.size() ||
        std::abs(ts_current.current() - _t_end) <
            std::numeric_limits<double>::epsilon())
    {
        return std::make_tuple(false, 0.0);
    }

    double dt = _dt_vector[ts_current.timeStepNumber()];
    if (ts_current.current() + dt > _t_end)
    {  // upper bound by t_end
        dt = _t_end - ts_current.current();
    }

    return std::make_tuple(true, dt);
}

}  // namespace NumLib
