/**
 * \author Norihiro Watanabe
 * \date   2012-08-03
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "FixedTimeStepping.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <limits>
#include <cassert>

#include <logog/include/logog.hpp>

#include "BaseLib/ConfigTree.h"
#include "BaseLib/Error.h"

namespace NumLib
{
FixedTimeStepping::FixedTimeStepping(double t0,
                                     double tn,
                                     const std::vector<double>& vec_all_dt)
    : _t_initial(t0),
      _t_end(computeEnd(t0, tn, vec_all_dt)),
      _dt_vector(vec_all_dt),
      _ts_prev(t0),
      _ts_current(t0)
{
}

FixedTimeStepping::FixedTimeStepping(double t0, double tn, double dt)
    : _t_initial(t0),
      _t_end(tn),
      _dt_vector(static_cast<std::size_t>(std::ceil((tn - t0) / dt)), dt),
      _ts_prev(t0),
      _ts_current(t0)
{
}

std::unique_ptr<ITimeStepAlgorithm> FixedTimeStepping::newInstance(
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

const TimeStep FixedTimeStepping::getTimeStep() const
{
    return _ts_current;
}

bool FixedTimeStepping::next()
{
    // check if last time step
    if (_ts_current.steps() == _dt_vector.size() ||
        std::abs(_ts_current.current() - _t_end) <
            std::numeric_limits<double>::epsilon())
        return false;

    // confirm current time and move to the next if accepted
    if (accepted())
        _ts_prev = _ts_current;

    // prepare the next time step info
    _ts_current = _ts_prev;
    double dt = _dt_vector[_ts_prev.steps()];
    if (_ts_prev.current() + dt > _t_end)  // upper bound by t_end
        dt = _t_end - _ts_prev.current();
    _ts_current += dt;

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

}  // NumLib
