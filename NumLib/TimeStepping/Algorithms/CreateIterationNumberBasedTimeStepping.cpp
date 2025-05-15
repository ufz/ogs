/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
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

IterationNumberBasedTimeSteppingParameters
parseIterationNumberBasedTimeStepping(BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{prj__time_loop__processes__process__time_stepping__type}
    config.checkConfigParameter("type", "IterationNumberBasedTimeStepping");

    //! \ogs_file_param{prj__time_loop__processes__process__time_stepping__IterationNumberBasedTimeStepping__t_initial}
    auto const t_initial = config.getConfigParameter<double>("t_initial");
    //! \ogs_file_param{prj__time_loop__processes__process__time_stepping__IterationNumberBasedTimeStepping__t_end}
    auto const t_end = config.getConfigParameter<double>("t_end");
    if (t_end < t_initial)
    {
        OGS_FATAL(
            "iteration number based timestepping: t_end({}) is smaller than "
            "t_initial({})",
            t_end,
            t_initial);
    }

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

    std::string const multiplier_interpolation_type_string =
        //! \ogs_file_param{prj__time_loop__processes__process__time_stepping__IterationNumberBasedTimeStepping__multiplier_interpolation_type}
        config.getConfigParameter<std::string>("multiplier_interpolation_type",
                                               "PiecewiseConstant");
    auto const multiplier_interpolation_type =
        (multiplier_interpolation_type_string == "PiecewiseLinear")
            ? MultiplyerInterpolationType::PiecewiseLinear
            : MultiplyerInterpolationType::PiecewiseConstant;

    return {t_initial,
            t_end,
            minimum_dt,
            maximum_dt,
            initial_dt,
            multiplier_interpolation_type,
            std::move(number_iterations),
            std::move(multiplier)};
}

/// Create a IterationNumberBasedTimeStepping time stepper from the given
/// configuration.
std::unique_ptr<TimeStepAlgorithm> createIterationNumberBasedTimeStepping(
    IterationNumberBasedTimeSteppingParameters&& parameters,
    std::vector<double> const& fixed_times_for_output)
{
    if (parameters.t_end < parameters.t_initial)
    {
        OGS_FATAL(
            "iteration number based timestepping: end time ({}) is smaller "
            "than initial time ({})",
            parameters.t_end,
            parameters.t_initial);
    }

    return std::make_unique<IterationNumberBasedTimeStepping>(
        parameters.t_initial, parameters.t_end, parameters.minimum_dt,
        parameters.maximum_dt, parameters.initial_dt,
        parameters.multiplier_interpolation_type,
        std::move(parameters.number_iterations),
        std::move(parameters.multiplier), fixed_times_for_output);
}

}  // namespace NumLib
