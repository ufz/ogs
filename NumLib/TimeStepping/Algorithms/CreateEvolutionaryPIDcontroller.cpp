/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  Created on June 26, 2017, 4:43 PM
 */

#include "CreateEvolutionaryPIDcontroller.h"

#include "BaseLib/Algorithm.h"
#include "BaseLib/ConfigTree.h"
#include "EvolutionaryPIDcontroller.h"
#include "TimeStepAlgorithm.h"

namespace NumLib
{
class TimeStepAlgorithm;

EvolutionaryPIDcontrollerParameters parseEvolutionaryPIDcontroller(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{prj__time_loop__processes__process__time_stepping__type}
    config.checkConfigParameter("type", "EvolutionaryPIDcontroller");

    //! \ogs_file_param{prj__time_loop__processes__process__time_stepping__EvolutionaryPIDcontroller__t_initial}
    auto const t0 = config.getConfigParameter<double>("t_initial");
    //! \ogs_file_param{prj__time_loop__processes__process__time_stepping__EvolutionaryPIDcontroller__t_end}
    auto const t_end = config.getConfigParameter<double>("t_end");

    //! \ogs_file_param{prj__time_loop__processes__process__time_stepping__EvolutionaryPIDcontroller__dt_guess}
    auto const h0 = config.getConfigParameter<double>("dt_guess");

    //! \ogs_file_param{prj__time_loop__processes__process__time_stepping__EvolutionaryPIDcontroller__dt_min}
    auto const h_min = config.getConfigParameter<double>("dt_min");
    //! \ogs_file_param{prj__time_loop__processes__process__time_stepping__EvolutionaryPIDcontroller__dt_max}
    auto const h_max = config.getConfigParameter<double>("dt_max");
    //! \ogs_file_param{prj__time_loop__processes__process__time_stepping__EvolutionaryPIDcontroller__rel_dt_min}
    auto const rel_h_min = config.getConfigParameter<double>("rel_dt_min");
    //! \ogs_file_param{prj__time_loop__processes__process__time_stepping__EvolutionaryPIDcontroller__rel_dt_max}
    auto const rel_h_max = config.getConfigParameter<double>("rel_dt_max");

    //! \ogs_file_param{prj__time_loop__processes__process__time_stepping__EvolutionaryPIDcontroller__tol}
    auto const tol = config.getConfigParameter<double>("tol");

    return {t0, t_end, h0, h_min, h_max, rel_h_min, rel_h_max, tol};
}

std::unique_ptr<TimeStepAlgorithm> createEvolutionaryPIDcontroller(
    EvolutionaryPIDcontrollerParameters const& config,
    std::vector<double> const& fixed_times_for_output)
{
    if (config.t_end < config.t0)
    {
        OGS_FATAL(
            "Evolutionary PID controller timestepping: end time ({}) is "
            "smaller than initial time ({})",
            config.t_end,
            config.t0);
    }

    return std::make_unique<EvolutionaryPIDcontroller>(config.t0,
                                                       config.t_end,
                                                       config.h0,
                                                       config.h_min,
                                                       config.h_max,
                                                       config.rel_h_min,
                                                       config.rel_h_max,
                                                       config.tol,
                                                       fixed_times_for_output);
}
}  // end of namespace NumLib
