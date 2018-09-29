/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  \file   CreateEvolutionaryPIDcontroller.cpp
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
std::unique_ptr<TimeStepAlgorithm> createEvolutionaryPIDcontroller(
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

    auto fixed_output_times =
        //! \ogs_file_param{prj__time_loop__processes__process__time_stepping__EvolutionaryPIDcontroller__fixed_output_times}
        config.getConfigParameter<std::vector<double>>("fixed_output_times",
                                                       std::vector<double>{});
    if (!fixed_output_times.empty())
    {
        // Remove possible duplicated elements and sort in descending order.
        BaseLib::makeVectorUnique(fixed_output_times, std::greater<double>());
    }

    //! \ogs_file_param{prj__time_loop__processes__process__time_stepping__EvolutionaryPIDcontroller__tol}
    auto const tol = config.getConfigParameter<double>("tol");

    return std::make_unique<EvolutionaryPIDcontroller>(
        t0, t_end, h0, h_min, h_max, rel_h_min, rel_h_max,
        std::move(fixed_output_times), tol);
}
}  // end of namespace NumLib
