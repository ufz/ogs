/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  Created on June 26, 2017, 5:03 PM
 */

#include "CreateFixedTimeStepping.h"

#include "BaseLib/ConfigTree.h"
#include "BaseLib/Error.h"
#include "FixedTimeStepping.h"
#include "TimeStepAlgorithm.h"

namespace NumLib
{
FixedTimeSteppingParameters parseFixedTimeStepping(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{prj__time_loop__processes__process__time_stepping__type}
    config.checkConfigParameter("type", "FixedTimeStepping");

    //! \ogs_file_param{prj__time_loop__processes__process__time_stepping__FixedTimeStepping__t_initial}
    auto const t_initial = config.getConfigParameter<double>("t_initial");
    //! \ogs_file_param{prj__time_loop__processes__process__time_stepping__FixedTimeStepping__t_end}
    auto const t_end = config.getConfigParameter<double>("t_end");

    //! \ogs_file_param{prj__time_loop__processes__process__time_stepping__FixedTimeStepping__n_steps}
    auto const n_steps = config.getConfigParameterOptional<int>("n_steps");
    if (n_steps.has_value())
    {
        if (t_end <= t_initial)
        {
            OGS_FATAL(
                "Creating linearly spaced time steps vector using "
                "FixedTimeStepping algorithm failed! "
                "User provided start value (t_initial) "
                "{} is not smaller then end value (t_end) {}.",
                t_initial, t_end);
        }

        if (*n_steps <= 0)
        {
            OGS_FATAL(
                "Requested number of time steps in time steps vector "
                "(n_steps) must be greater then 0. "
                "{} time steps were requested",
                *n_steps);
        }
        // Create the RepeatDtPair
        double const t_step =
            (t_end - t_initial) / static_cast<double>(*n_steps);
        std::vector const repeat_pairs = {
            RepeatDtPair{static_cast<std::size_t>(*n_steps), t_step}};
        return {t_initial, t_end, repeat_pairs};
    }

    //! \ogs_file_param{prj__time_loop__processes__process__time_stepping__FixedTimeStepping__timesteps}
    auto const delta_ts_config = config.getConfigSubtree("timesteps");

    // TODO: consider adding call "listNonEmpty" to config tree
    //! \ogs_file_param{prj__time_loop__processes__process__time_stepping__FixedTimeStepping__timesteps__pair}
    auto const range = delta_ts_config.getConfigSubtreeList("pair");
    if (range.begin() == range.end())
    {
        OGS_FATAL("no timesteps have been given");
    }

    std::vector<RepeatDtPair> repeat_dt_pairs;
    for (auto const pair : range)
    {
        repeat_dt_pairs.emplace_back(
            //! \ogs_file_param{prj__time_loop__processes__process__time_stepping__FixedTimeStepping__timesteps__pair__repeat}
            pair.getConfigParameter<std::size_t>("repeat"),
            //! \ogs_file_param{prj__time_loop__processes__process__time_stepping__FixedTimeStepping__timesteps__pair__delta_t}
            pair.getConfigParameter<double>("delta_t"));
    }

    return {t_initial, t_end, repeat_dt_pairs};
}

std::unique_ptr<TimeStepAlgorithm> createFixedTimeStepping(
    FixedTimeSteppingParameters const& parameters,
    std::vector<double> const& fixed_times_for_output)
{
    if (parameters.t_end < parameters.t_initial)
    {
        OGS_FATAL(
            "fixed timestepping: end time ({}) is smaller than initial time "
            "({})",
            parameters.t_end,
            parameters.t_initial);
    }

    if (!FixedTimeStepping::areRepeatDtPairsValid(parameters.repeat_dt_pairs))
    {
        OGS_FATAL(
            "CreateFixedTimeStepping: invalid specification of (repeat, "
            "delta_t) pairs");
    }

    return std::make_unique<FixedTimeStepping>(parameters.t_initial,
                                               parameters.t_end,
                                               parameters.repeat_dt_pairs,
                                               fixed_times_for_output);
}
}  // end of namespace NumLib
