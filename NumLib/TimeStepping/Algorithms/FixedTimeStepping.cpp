/**
 * \author Norihiro Watanabe
 * \date   2012-08-03
 *
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
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

#include "BaseLib/ConfigTree.h"
#include "logog/include/logog.hpp"

namespace NumLib
{

FixedTimeStepping::FixedTimeStepping(double t0, double tn, const std::vector<double> &vec_all_dt)
: _t_initial(t0), _t_end(computeEnd(t0, tn, vec_all_dt)), _dt_vector(vec_all_dt), _ts_prev(t0), _ts_current(t0)
{}

FixedTimeStepping::FixedTimeStepping(double t0, double tn, double dt)
: _t_initial(t0), _t_end(tn), _dt_vector(static_cast<std::size_t>(std::ceil((tn-t0)/dt)), dt), _ts_prev(t0), _ts_current(t0)
{}

FixedTimeStepping*
FixedTimeStepping::newInstance(BaseLib::ConfigTree const& config)
{
    assert(config.get<std::string>("type") == "FixedTimeStepping");

    auto const t_initial = config.get_optional<double>("t_initial");
    auto const t_end     = config.get_optional<double>("t_end");
    auto const dt        = config.get_optional<double>("dt");

    if (!t_initial)
    {
        ERR("could not find required parameter t_initial.");
        return nullptr;
    }
    if (!t_end)
    {
        ERR("could not find required parameter t_end.");
        return nullptr;
    }
    if (!dt)
    {
        ERR("could not find required parameter dt.");
        return nullptr;
    }

    return new FixedTimeStepping(*t_initial, *t_end, *dt);
}

const TimeStep FixedTimeStepping::getTimeStep() const
{
    return _ts_current;
}

bool FixedTimeStepping::next()
{
    // check if last time step
    if (_ts_current.steps() == _dt_vector.size()
        || std::abs(_ts_current.current()-_t_end) < std::numeric_limits<double>::epsilon())
        return false;

    // confirm current time and move to the next if accepted
    if (accepted())
        _ts_prev = _ts_current;

    // prepare the next time step info
    _ts_current = _ts_prev;
    double dt = _dt_vector[_ts_prev.steps()];
    if (_ts_prev.current() + dt > _t_end) // upper bound by t_end
        dt = _t_end - _ts_prev.current();
    _ts_current += dt;

    return true;
}

double FixedTimeStepping::computeEnd(double t_initial, double t_end, const std::vector<double> &dt_vector)
{
    double t_sum = t_initial + std::accumulate(dt_vector.begin(), dt_vector.end(), 0.);
    return std::min(t_end, t_sum);
}

} //NumLib
