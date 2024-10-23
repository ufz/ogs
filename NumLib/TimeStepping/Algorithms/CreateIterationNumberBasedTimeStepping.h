/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <memory>
#include <vector>

namespace BaseLib
{
class ConfigTree;
}

namespace NumLib
{
class TimeStepAlgorithm;

struct IterationNumberBasedTimeSteppingParameters final
{
    double t_initial;
    double t_end;
    double minimum_dt;
    double maximum_dt;
    double initial_dt;
    std::vector<int> number_iterations;
    std::vector<double> multiplier;
};

IterationNumberBasedTimeSteppingParameters
parseIterationNumberBasedTimeStepping(BaseLib::ConfigTree const& config);

/// Create a IterationNumberBasedTimeStepping time stepper from the given
/// configuration.
std::unique_ptr<TimeStepAlgorithm> createIterationNumberBasedTimeStepping(
    IterationNumberBasedTimeSteppingParameters&& parameters,
    std::vector<double> const& fixed_times_for_output);
}  // namespace NumLib
