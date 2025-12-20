// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <memory>
#include <vector>

#include "MultiplyerInterpolationType.h"

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
    MultiplyerInterpolationType multiplier_interpolation_type;
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
