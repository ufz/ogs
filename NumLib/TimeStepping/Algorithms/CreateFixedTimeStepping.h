/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  Created on June 26, 2017, 5:02 PM
 */

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
