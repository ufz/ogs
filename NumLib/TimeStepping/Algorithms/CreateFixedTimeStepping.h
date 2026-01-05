// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <memory>
#include <vector>

#include "FixedTimeStepping.h"

namespace BaseLib
{
class ConfigTree;
}

namespace NumLib
{
class TimeStepAlgorithm;

struct FixedTimeSteppingParameters final
{
    double t_initial;
    double t_end;
    std::vector<RepeatDtPair> repeat_dt_pairs;
};

/// Create a FixedTimeStepping time stepper from the given
/// configuration
FixedTimeSteppingParameters parseFixedTimeStepping(
    BaseLib::ConfigTree const& config);

std::unique_ptr<TimeStepAlgorithm> createFixedTimeStepping(
    FixedTimeSteppingParameters const& parameters,
    std::vector<double> const& fixed_times_for_output);

}  // end of namespace NumLib
