/**
 * \file
 * \author Norihiro Watanabe
 * \date   2012-08-03
 *
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "FixedTimeStepping.h"

#include <algorithm>
#include <cassert>
#include <limits>
#include <numeric>

namespace
{
/// Returns sum of the newly added time increments.
double addTimeIncrement(std::vector<double>& delta_ts, std::size_t const repeat,
                        double const delta_t, double const t_curr)
{
    auto const new_size = delta_ts.size() + repeat;
    try
    {
        delta_ts.resize(new_size, delta_t);
    }
    catch (std::exception const& e)
    {
        OGS_FATAL(
            "Resize of the time steps vector failed for the requested new size "
            "{:d}. Probably there is not enough memory ({:g} GiB "
            "requested).\nThrown exception: {:s}",
            new_size,
            new_size * sizeof(double) / 1024. / 1024. / 1024.,
            e.what());
    }

    // Multiplying dt * repeat is not the same as in the current
    // implementation of the time loop, where the dt's are added.
    // Therefore the sum of all dt is taken here.
    return std::accumulate(end(delta_ts) - repeat, end(delta_ts), t_curr);
}
}  // namespace

namespace NumLib
{
std::size_t findDeltatInterval(double const t_initial,
                               std::vector<double> const& delta_ts,
                               double const fixed_output_time)
{
    if (fixed_output_time < t_initial)
    {
        return std::numeric_limits<std::size_t>::max();
    }

    auto timestepper_time = t_initial;
    for (std::size_t k = 0; k < delta_ts.size(); ++k)
    {
        if (timestepper_time <= fixed_output_time &&
            fixed_output_time < timestepper_time + delta_ts[k])
        {
            return k;
        }
        timestepper_time += delta_ts[k];
    }

    return std::numeric_limits<std::size_t>::max();
}

void incorporateFixedTimesForOutput(
    double const t_initial, double const t_end, std::vector<double>& delta_ts,
    std::vector<double> const& fixed_times_for_output)
{
    if (fixed_times_for_output.empty())
    {
        return;
    }

    if (auto lower_bound =
            std::lower_bound(begin(fixed_times_for_output),
                             end(fixed_times_for_output), t_initial);
        lower_bound != begin(fixed_times_for_output))
    {
        WARN(
            "Request for output at times {}, but the simulation's start time "
            "is {}. Output will be skipped.",
            fmt::join(begin(fixed_times_for_output), lower_bound, ", "),
            t_initial);
    }

    if (auto upper_bound = std::upper_bound(begin(fixed_times_for_output),
                                            end(fixed_times_for_output), t_end);
        upper_bound != end(fixed_times_for_output))
    {
        WARN(
            "Request for output at times {}, but simulation's end time is {}. "
            "Output will be skipped.",
            fmt::join(upper_bound, end(fixed_times_for_output), ", "),
            t_end);
    }

    if (delta_ts.empty())
    {
        WARN("No timesteps specified.");
        return;
    }

    // incorporate fixed output times into dts vector
    for (auto const fixed_time_for_output : fixed_times_for_output)
    {
        auto const interval_number =
            findDeltatInterval(t_initial, delta_ts, fixed_time_for_output);
        if (interval_number == std::numeric_limits<std::size_t>::max())
        {
            WARN("Did not find interval for fixed output time {}",
                 fixed_time_for_output);
            continue;
        }
        auto const lower_bound = std::accumulate(
            begin(delta_ts), begin(delta_ts) + interval_number, t_initial);
        auto const upper_bound = lower_bound + delta_ts[interval_number];
        if (fixed_time_for_output - lower_bound <=
            std::numeric_limits<double>::epsilon())
        {
            continue;
        }
        if (upper_bound - fixed_time_for_output <=
            std::numeric_limits<double>::epsilon())
        {
            continue;
        }
        delta_ts[interval_number] = fixed_time_for_output - lower_bound;

        delta_ts.insert(delta_ts.begin() + interval_number + 1,
                        upper_bound - fixed_time_for_output);
    }
}

FixedTimeStepping::FixedTimeStepping(
    double t0, double tn,
    std::vector<std::pair<std::size_t, double>> const& repeat_dt_pairs,
    std::vector<double> const& fixed_times_for_output)
    : TimeStepAlgorithm(t0, tn)
{
    double t_curr = _t_initial;

    for (auto const& [repeat, delta_t] : repeat_dt_pairs)
    {
        if (repeat == 0)
        {
            OGS_FATAL("<repeat> is zero.");
        }
        if (delta_t <= 0.0)
        {
            OGS_FATAL("timestep <delta_t> is <= 0.0.");
        }

        if (t_curr <= _t_end)
        {
            t_curr = addTimeIncrement(_dt_vector, repeat, delta_t, t_curr);
        }
    }

    // append last delta_t until t_end is reached
    if (t_curr <= _t_end)
    {
        auto const delta_t = repeat_dt_pairs.back().second;
        auto const repeat =
            static_cast<std::size_t>(std::ceil((_t_end - t_curr) / delta_t));
        addTimeIncrement(_dt_vector, repeat, delta_t, t_curr);
    }

    incorporateFixedTimesForOutput(_t_initial, _t_end, _dt_vector,
                                   fixed_times_for_output);
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
            "size {}. Probably there is not enough memory ({:g} GiB "
            "requested).\n"
            "Thrown exception: {}",
            new_size, new_size * sizeof(double) / 1024. / 1024. / 1024.,
            e.what());
    }
    catch (std::bad_alloc const& e)
    {
        OGS_FATAL(
            "Allocation of the time steps vector failed for the requested "
            "size {}. Probably there is not enough memory ({:g} GiB "
            "requested).\n"
            "Thrown exception: {}",
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
        std::abs(ts_current.current() - end()) <
            std::numeric_limits<double>::epsilon())
    {
        return std::make_tuple(false, 0.0);
    }

    double dt = _dt_vector[ts_current.timeStepNumber()];
    if (ts_current.current() + dt > end())
    {  // upper bound by t_end
        dt = end() - ts_current.current();
    }

    return std::make_tuple(true, dt);
}

}  // namespace NumLib
