/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  \file   EvolutionaryPIDcontroller.cpp
 *  Created on March 31, 2017, 4:13 PM
 */

#include <algorithm>
#include <functional>
#include <limits>
#include <vector>
#include <logog/include/logog.hpp>

#include "BaseLib/makeVectorUnique.h"

#include "EvolutionaryPIDcontroller.h"

namespace NumLib
{
bool EvolutionaryPIDcontroller::next(const double solution_error)
{
    const bool is_previous_step_accepted = _is_accepted;

    const double e_n = solution_error;
    const double zero_threshlod = std::numeric_limits<double>::epsilon();
    // step rejected.
    if (e_n > _tol)  // e_n < TOL
    {
        _is_accepted = false;

        double h_new = (e_n > zero_threshlod) ? _ts_current.dt() * _tol / e_n
                                              : 0.5 * _ts_current.dt();

        h_new = limitStepSize(h_new, is_previous_step_accepted);
        h_new = checkSpecificTimeReached(h_new);

        _ts_current = _ts_prev;
        _ts_current += h_new;

        WARN(
            "This step is rejected due to the relative change from the"
            " solution of the previous\n"
            "\t time step to the current solution exceeds the given tolerance"
            " of %g.\n"
            "\t This time step will be repeated with a new time step size of"
            " %g\n"
            "\t or the simulation will be halted.",
            _tol, h_new);

        return false;
    }

    // step accepted.
    _is_accepted = true;

    if (_ts_current.steps() == 0)
    {
        _ts_prev = _ts_current;
        _ts_current += _h0;
        _e_n_minus1 = e_n;

        _dt_vector.push_back(_h0);
    }
    else
    {
        const double h_n = _ts_current.dt();
        double h_new = h_n;

        if (e_n > zero_threshlod)
        {
            if (_e_n_minus1 > zero_threshlod)
            {
                if (_e_n_minus2 > zero_threshlod)
                {
                    h_new = std::pow(_e_n_minus1 / e_n, _kP) *
                            std::pow(_tol / e_n, _kI) *
                            std::pow(
                                _e_n_minus1 * _e_n_minus1 / (e_n * _e_n_minus2),
                                _kD) *
                            h_n;
                }
                else
                {
                    h_new = std::pow(_e_n_minus1 / e_n, _kP) *
                            std::pow(_tol / e_n, _kI) * h_n;
                }
            }
            else
            {
                h_new = std::pow(_tol / e_n, _kI) * h_n;
            }
        }

        h_new = limitStepSize(h_new, is_previous_step_accepted);
        h_new = checkSpecificTimeReached(h_new);
        _dt_vector.push_back(h_new);

        _ts_prev = _ts_current;
        _ts_current += h_new;

        _e_n_minus2 = _e_n_minus1;
        _e_n_minus1 = e_n;
    }

    return true;
}

double EvolutionaryPIDcontroller::limitStepSize(
    const double h_new, const bool previous_step_accepted) const
{
    const double h_n = _ts_current.dt();
    // Forced the computed time step size in the given range
    // (see the formulas in the documentation of the class)
    const double h_in_range = std::max(_h_min, std::min(h_new, _h_max));
    // Limit the step size change ratio.
    double limited_h =
        std::max(_rel_h_min * h_n, std::min(h_in_range, _rel_h_max * h_n));

    if (!previous_step_accepted)
    {
        // If the last time step was rejected and the new predicted time step
        // size is identical to that of the previous rejected step, the new
        // step size is then reduced by half.
        if (std::fabs(limited_h - _ts_current.dt()) <
            std::numeric_limits<double>::min())
            limited_h = std::max(_h_min, 0.5 * limited_h);

        // If the last time step was rejected and the new predicted time step
        // size is larger than the step size of the rejected step, the new step
        // size takes the half of the size of the rejected step. This could
        // happen when a time step is rejected due to a diverged non-linear
        // solver. In such case, this algorithm may give a large time step size
        // by using the diverged solution.
        if (limited_h > _ts_current.dt())
            limited_h = std::max(_h_min, 0.5 * _ts_current.dt());
    }
    return limited_h;
}

double EvolutionaryPIDcontroller::checkSpecificTimeReached(const double h_new)
{
    if (_specific_times.empty())
        return h_new;

    const double specific_time = _specific_times.back();
    const double zero_threshold = std::numeric_limits<double>::epsilon();
    if ((specific_time > _ts_current.current()) &&
        (_ts_current.current() + h_new - specific_time > zero_threshold))
    {
        _specific_times.pop_back();
        return specific_time - _ts_current.current();
    }

    return h_new;
}

void EvolutionaryPIDcontroller::addSpecificTimes(
    std::vector<double> const& extra_specific_times)
{
    _specific_times.insert(_specific_times.end(), extra_specific_times.begin(),
                           extra_specific_times.end());

    // Remove possible duplicated elements and sort in descending order.
    BaseLib::makeVectorUnique(_specific_times, std::greater<double>());
}

}  // end of namespace NumLib
