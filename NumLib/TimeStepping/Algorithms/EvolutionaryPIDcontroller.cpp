/**
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  \file
 *  Created on March 31, 2017, 4:13 PM
 */

#include "EvolutionaryPIDcontroller.h"

#include <functional>
#include <limits>
#include <vector>

#include "BaseLib/Algorithm.h"
#include "BaseLib/Logging.h"

namespace NumLib
{
std::tuple<bool, double> EvolutionaryPIDcontroller::next(
    double const solution_error, int const /*number_iterations*/)
{
    const bool is_previous_step_accepted = _is_accepted;

    const double e_n = solution_error;
    const double zero_threshlod = std::numeric_limits<double>::epsilon();
    // step rejected.
    if (e_n > _tol)
    {
        _is_accepted = false;

        double h_new = (e_n > zero_threshlod) ? _ts_current.dt() * _tol / e_n
                                              : 0.5 * _ts_current.dt();

        h_new = limitStepSize(h_new, is_previous_step_accepted);

        WARN(
            "This step is rejected due to the relative change from the"
            " solution of the previous\n"
            "\t time step to the current solution exceeds the given tolerance"
            " of {:g}.\n"
            "\t This time step will be repeated with a new time step size of"
            " {:g}\n"
            "\t or the simulation will be halted.",
            _tol, h_new);

        return std::make_tuple(false, h_new);
    }

    // step accepted.
    _is_accepted = true;

    if (_ts_current.timeStepNumber() == 0)
    {
        _e_n_minus1 = e_n;

        return std::make_tuple(true, _h0);
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

        _e_n_minus2 = _e_n_minus1;
        _e_n_minus1 = e_n;

        return std::make_tuple(true, h_new);
    }

    return {};
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
        if (std::abs(limited_h - _ts_current.dt()) <
            std::numeric_limits<double>::min())
        {
            limited_h = std::max(_h_min, 0.5 * limited_h);
        }

        // If the last time step was rejected and the new predicted time step
        // size is larger than the step size of the rejected step, the new step
        // size takes the half of the size of the rejected step. This could
        // happen when a time step is rejected due to a diverged non-linear
        // solver. In such case, this algorithm may give a large time step size
        // by using the diverged solution.
        if (limited_h > _ts_current.dt())
        {
            limited_h = std::max(_h_min, 0.5 * _ts_current.dt());
        }
    }
    return limited_h;
}

bool EvolutionaryPIDcontroller::canReduceTimestepSize() const
{
    // If current and previous dt are both at minimum dt, then cannot reduce
    // further.
    return !(_ts_current.dt() == _h_min && _ts_prev.dt() == _h_min);
}
}  // namespace NumLib
