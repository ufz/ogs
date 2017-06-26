/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  \file   CreateFixedTimeStepping.cpp
 *  Created on June 26, 2017, 5:03 PM
 */

#include "CreateFixedTimeStepping.h"
#include <string>

#include "BaseLib/ConfigTree.h"
#include "BaseLib/Error.h"

#include "FixedTimeStepping.h"
#include "TimeStepAlgorithm.h"

namespace NumLib
{
class TimeStepAlgorithm;
std::unique_ptr<TimeStepAlgorithm> createFixedTimeStepping(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{prj__time_loop__time_stepping__type}
    config.checkConfigParameter("type", "FixedTimeStepping");

    //! \ogs_file_param{prj__time_loop__time_stepping__FixedTimeStepping__t_initial}
    auto const t_initial = config.getConfigParameter<double>("t_initial");
    //! \ogs_file_param{prj__time_loop__time_stepping__FixedTimeStepping__t_end}
    auto const t_end = config.getConfigParameter<double>("t_end");
    //! \ogs_file_param{prj__time_loop__time_stepping__FixedTimeStepping__timesteps}
    auto const delta_ts = config.getConfigSubtree("timesteps");

    std::vector<double> timesteps;
    double t_curr = t_initial;
    double delta_t = 0.0;

    // TODO: consider adding call "listNonEmpty" to config tree
    //! \ogs_file_param{prj__time_loop__time_stepping__FixedTimeStepping__timesteps__pair}
    auto const range = delta_ts.getConfigSubtreeList("pair");
    if (range.begin() == range.end())
    {
        OGS_FATAL("no timesteps have been given");
    }
    for (auto const pair : range)
    {
        //! \ogs_file_param{prj__time_loop__time_stepping__FixedTimeStepping__timesteps__pair__repeat}
        auto const repeat = pair.getConfigParameter<std::size_t>("repeat");
        //! \ogs_file_param{prj__time_loop__time_stepping__FixedTimeStepping__timesteps__pair__delta_t}
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
            timesteps.resize(timesteps.size() + repeat, delta_t);

            t_curr += repeat * delta_t;
        }
    }

    // append last delta_t until t_end is reached
    if (t_curr <= t_end)
    {
        auto const repeat =
            static_cast<std::size_t>(std::ceil((t_end - t_curr) / delta_t));
        timesteps.resize(timesteps.size() + repeat, delta_t);
    }

    return std::make_unique<FixedTimeStepping>(t_initial, t_end, timesteps);
}
}  // end of namespace NumLib
