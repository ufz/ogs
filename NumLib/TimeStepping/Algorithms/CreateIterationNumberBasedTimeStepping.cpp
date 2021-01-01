/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "CreateIterationNumberBasedTimeStepping.h"

#include <string>

#include "BaseLib/Algorithm.h"
#include "BaseLib/ConfigTree.h"
#include "IterationNumberBasedTimeStepping.h"
#include "TimeStepAlgorithm.h"

namespace NumLib
{
class TimeStepAlgorithm;
std::unique_ptr<TimeStepAlgorithm> createIterationNumberBasedTimeStepping(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{prj__time_loop__processes__process__time_stepping__type}
    config.checkConfigParameter("type", "IterationNumberBasedTimeStepping");

    //! \ogs_file_param{prj__time_loop__processes__process__time_stepping__IterationNumberBasedTimeStepping__t_initial}
    auto const t_initial = config.getConfigParameter<double>("t_initial");
    //! \ogs_file_param{prj__time_loop__processes__process__time_stepping__IterationNumberBasedTimeStepping__t_end}
    auto const t_end = config.getConfigParameter<double>("t_end");
    //! \ogs_file_param{prj__time_loop__processes__process__time_stepping__IterationNumberBasedTimeStepping__initial_dt}
    auto const initial_dt = config.getConfigParameter<double>("initial_dt");
    //! \ogs_file_param{prj__time_loop__processes__process__time_stepping__IterationNumberBasedTimeStepping__minimum_dt}
    auto const minimum_dt = config.getConfigParameter<double>("minimum_dt");
    //! \ogs_file_param{prj__time_loop__processes__process__time_stepping__IterationNumberBasedTimeStepping__maximum_dt}
    auto const maximum_dt = config.getConfigParameter<double>("maximum_dt");

    auto number_iterations =
        //! \ogs_file_param{prj__time_loop__processes__process__time_stepping__IterationNumberBasedTimeStepping__number_iterations}
        config.getConfigParameter<std::vector<int>>("number_iterations");
    auto multiplier =
        //! \ogs_file_param{prj__time_loop__processes__process__time_stepping__IterationNumberBasedTimeStepping__multiplier}
        config.getConfigParameter<std::vector<double>>("multiplier");

    return std::make_unique<IterationNumberBasedTimeStepping>(
        t_initial, t_end, minimum_dt, maximum_dt, initial_dt,
        std::move(number_iterations), std::move(multiplier));
}
}  // namespace NumLib
