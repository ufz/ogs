/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  Created on May 2, 2017, 12:18 PM
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
std::unique_ptr<TimeStepAlgorithm> createTimeStepper(
    BaseLib::ConfigTree const& config,
    std::vector<double> const& fixed_times_for_output);
};  // namespace NumLib
