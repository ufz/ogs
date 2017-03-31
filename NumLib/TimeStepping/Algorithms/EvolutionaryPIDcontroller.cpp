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
#include <string>
#include <vector>

#include "EvolutionaryPIDcontroller.h"

#include "BaseLib/ConfigTree.h"
#include "ITimeStepAlgorithm.h"

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

        const double h_new = (e_n > zero_threshlod && _tol / e_n < 1.0)
                                 ? _ts_prev.dt() * _tol / e_n
                                 : 0.5 * _ts_prev.dt();

        _ts_current = _ts_prev;
        _ts_current += h_new;
        return false;
    }

    // step accepted.
    _is_accepted = true;
    if (_specific_times.size() > 0)
        _specific_times.pop_back();

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

            h_new = limitStepSize(h_new, h_n);
        }

        _dt_vector.push_back(h_new);

        _ts_prev = _ts_current;
        _ts_current += h_new;

        _e_n_minus2 = _e_n_minus1;
        _e_n_minus1 = e_n;
    }

    return true;
}

/// Create an EvolutionaryPIDcontroller time stepper from the given
/// configuration
std::unique_ptr<ITimeStepAlgorithm> createEvolutionaryPIDcontroller(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{prj__time_loop__time_stepping__type}
    config.checkConfigParameter("type", "EvolutionaryPIDcontroller");

    //! \ogs_file_param{prj__time_loop__time_stepping__EvolutionaryPIDcontroller__t_initial}
    auto const t0 = config.getConfigParameter<double>("t_initial");
    //! \ogs_file_param{prj__time_loop__time_stepping__EvolutionaryPIDcontroller__t_end}
    auto const t_end = config.getConfigParameter<double>("t_end");
    //! \ogs_file_param{prj__time_loop__time_stepping__EvolutionaryPIDcontroller__dt_guess}
    auto const h0 = config.getConfigParameter<double>("dt_guess");

    //! \ogs_file_param{prj__time_loop__time_stepping__EvolutionaryPIDcontroller__dt_min}
    auto const h_min = config.getConfigParameter<double>("dt_min");
    //! \ogs_file_param{prj__time_loop__time_stepping__EvolutionaryPIDcontroller__dt_max}
    auto const h_max = config.getConfigParameter<double>("dt_max");
    //! \ogs_file_param{prj__time_loop__time_stepping__EvolutionaryPIDcontroller__rel_dt_min}
    auto const rel_h_min = config.getConfigParameter<double>("rel_dt_min");
    //! \ogs_file_param{prj__time_loop__time_stepping__EvolutionaryPIDcontroller__rel_dt_max}
    auto const rel_h_max = config.getConfigParameter<double>("rel_dt_max");

    auto specific_times_opt =
        //! \ogs_file_param{prj__time_loop__time_stepping__EvolutionaryPIDcontroller__specific_times}
        config.getConfigParameterOptional<std::vector<double>>(
            "specific_times");
    std::vector<double> dummy_vector;
    std::vector<double>& specific_times =
        (specific_times_opt) ? *specific_times_opt : dummy_vector;
    if (specific_times.size() > 0)
    {
        // Sort in descending order.
        std::sort(specific_times.begin(),
                  specific_times.end(),
                  std::greater<double>());
        // Remove possible duplicated elements.
        auto last = std::unique(specific_times.begin(), specific_times.end());
        specific_times.erase(last, specific_times.end());
    }

    //! \ogs_file_param{prj__time_loop__time_stepping__EvolutionaryPIDcontroller__tol}
    auto const tol = config.getConfigParameter<double>("tol");
    auto const norm_type_opt =
        //! \ogs_file_param{prj__time_loop__time_stepping__EvolutionaryPIDcontroller__norm_type}
        config.getConfigParameterOptional<std::string>("norm_type");
    const MathLib::VecNormType norm_type =
        (norm_type_opt) ? MathLib::convertStringToVecNormType(*norm_type_opt)
                        : MathLib::VecNormType::NORM2;

    return std::unique_ptr<ITimeStepAlgorithm>(
        new EvolutionaryPIDcontroller(t0,
                                      t_end,
                                      h0,
                                      h_min,
                                      h_max,
                                      rel_h_min,
                                      rel_h_max,
                                      std::move(specific_times),
                                      tol,
                                      norm_type));
}

}  // end of namespace NumLib
