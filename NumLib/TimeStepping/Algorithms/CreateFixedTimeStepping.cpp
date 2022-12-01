/**
 * \file
 * \copyright
 * Copyright (c) 2012-2023, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  Created on June 26, 2017, 5:03 PM
 */

#include "CreateFixedTimeStepping.h"

#include <fmt/ranges.h>

#include <algorithm>
#include <numeric>
#include <string>

#include "BaseLib/Algorithm.h"
#include "BaseLib/ConfigTree.h"
#include "BaseLib/Error.h"
#include "FixedTimeStepping.h"
#include "TimeStepAlgorithm.h"

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
    if (fixed_output_time == timestepper_time + delta_ts.back())
    {
        return std::numeric_limits<std::size_t>::max();
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
        auto const upper_bound = std::accumulate(
            begin(delta_ts), begin(delta_ts) + interval_number + 1, t_initial);
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

class TimeStepAlgorithm;
std::unique_ptr<TimeStepAlgorithm> createFixedTimeStepping(
    BaseLib::ConfigTree const& config,
    std::vector<double> const& fixed_times_for_output)
{
    //! \ogs_file_param{prj__time_loop__processes__process__time_stepping__type}
    config.checkConfigParameter("type", "FixedTimeStepping");

    //! \ogs_file_param{prj__time_loop__processes__process__time_stepping__FixedTimeStepping__t_initial}
    auto const t_initial = config.getConfigParameter<double>("t_initial");
    //! \ogs_file_param{prj__time_loop__processes__process__time_stepping__FixedTimeStepping__t_end}
    auto const t_end = config.getConfigParameter<double>("t_end");
    //! \ogs_file_param{prj__time_loop__processes__process__time_stepping__FixedTimeStepping__timesteps}
    auto const delta_ts_config = config.getConfigSubtree("timesteps");

    std::vector<double> delta_ts;
    double t_curr = t_initial;
    double delta_t = 0.0;

    // TODO: consider adding call "listNonEmpty" to config tree
    //! \ogs_file_param{prj__time_loop__processes__process__time_stepping__FixedTimeStepping__timesteps__pair}
    auto const range = delta_ts_config.getConfigSubtreeList("pair");
    if (range.begin() == range.end())
    {
        OGS_FATAL("no timesteps have been given");
    }
    for (auto const pair : range)
    {
        //! \ogs_file_param{prj__time_loop__processes__process__time_stepping__FixedTimeStepping__timesteps__pair__repeat}
        auto const repeat = pair.getConfigParameter<std::size_t>("repeat");
        //! \ogs_file_param{prj__time_loop__processes__process__time_stepping__FixedTimeStepping__timesteps__pair__delta_t}
        delta_t = pair.getConfigParameter<double>("delta_t");

        if (repeat == 0)
        {
            OGS_FATAL("<repeat> is zero.");
        }
        if (delta_t <= 0.0)
        {
            OGS_FATAL("timestep <delta_t> is <= 0.0.");
        }

        if (t_curr <= t_end)
        {
            auto const new_size = delta_ts.size() + repeat;
            try
            {
                delta_ts.resize(new_size, delta_t);
            }
            catch (std::length_error const& e)
            {
                OGS_FATAL(
                    "Resize of the time steps vector failed for the requested "
                    "new size {:d}. Probably there is not enough memory ({:g} "
                    "GiB requested).\n"
                    "Thrown exception: {:s}",
                    new_size,
                    new_size * sizeof(double) / 1024. / 1024. / 1024.,
                    e.what());
            }
            catch (std::bad_alloc const& e)
            {
                OGS_FATAL(
                    "Resize of the time steps vector failed for the requested "
                    "new size {:d}. Probably there is not enough memory ({:g} "
                    "GiB requested).\n"
                    "Thrown exception: {:s}",
                    new_size, new_size * sizeof(double) / 1024. / 1024. / 1024.,
                    e.what());
            }

            t_curr += repeat * delta_t;
        }
    }

    // append last delta_t until t_end is reached
    if (t_curr <= t_end)
    {
        auto const repeat =
            static_cast<std::size_t>(std::ceil((t_end - t_curr) / delta_t));
        auto const new_size = delta_ts.size() + repeat;
        try
        {
            delta_ts.resize(new_size, delta_t);
        }
        catch (std::length_error const& e)
        {
            OGS_FATAL(
                "Resize of the time steps vector failed for the requested new "
                "size {:d}. Probably there is not enough memory ({:g} GiB "
                "requested).\n"
                "Thrown exception: {:s}",
                new_size,
                new_size * sizeof(double) / 1024. / 1024. / 1024.,
                e.what());
        }
        catch (std::bad_alloc const& e)
        {
            OGS_FATAL(
                "Resize of the time steps vector failed for the requested new "
                "size {:d}. Probably there is not enough memory ({:g} GiB "
                "requested).\n"
                "Thrown exception: {:s}",
                new_size, new_size * sizeof(double) / 1024. / 1024. / 1024.,
                e.what());
        }
    }

    incorporateFixedTimesForOutput(t_initial, t_end, delta_ts,
                                   fixed_times_for_output);
    return std::make_unique<FixedTimeStepping>(t_initial, t_end, delta_ts);
}
}  // end of namespace NumLib
