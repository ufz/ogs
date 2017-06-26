/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
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

#include "EvolutionaryPIDcontroller.h"

namespace NumLib
{
bool EvolutionaryPIDcontroller::next(const double solution_error)
{
    const double e_n = solution_error;
    const double zero_threshlod = std::numeric_limits<double>::epsilon();
    // step rejected.
    if (e_n > _tol)  // e_n < TOL
    {
        _is_accepted = false;

        const double h_new = (e_n > zero_threshlod)
                                 ? _ts_current.dt() * _tol / e_n
                                 : 0.5 * _ts_current.dt();

        _ts_current = _ts_prev;
        _ts_current += h_new;
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
                    h_new = std::pow(_e_n_minus1 / e_n, _kP) *
                            std::pow(_tol / e_n, _kI) *
                            std::pow(
                                _e_n_minus1 * _e_n_minus1 / (e_n * _e_n_minus2),
                                _kD) *
                            h_n;
                else
                    h_new = std::pow(_e_n_minus1 / e_n, _kP) *
                            std::pow(_tol / e_n, _kI) * h_n;
            }
            else
            {
                h_new = std::pow(_tol / e_n, _kI) * h_n;
            }
        }

        h_new = limitStepSize(h_new, h_n);
        const double checked_h_new = checkSpecificTimeReached(h_new);
        _dt_vector.push_back(checked_h_new);

        _ts_prev = _ts_current;
        _ts_current += checked_h_new;

        _e_n_minus2 = _e_n_minus1;
        _e_n_minus1 = e_n;
    }

    return true;
}

double EvolutionaryPIDcontroller::checkSpecificTimeReached(const double h_new)
{
    if (_specific_times.empty())
        return h_new;

    const double specific_time = _specific_times.back();
    const double zero_threshlod = std::numeric_limits<double>::epsilon();
    if ((specific_time > _ts_current.current()) &&
        (_ts_current.current() + h_new - specific_time > zero_threshlod))
    {
        _specific_times.pop_back();
        return specific_time - _ts_current.current();
    }

    return h_new;
}

}  // end of namespace NumLib
