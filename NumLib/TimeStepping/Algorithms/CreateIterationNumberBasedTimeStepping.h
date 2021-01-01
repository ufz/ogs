/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <memory>

namespace BaseLib
{
class ConfigTree;
}

namespace NumLib
{
class TimeStepAlgorithm;
}

namespace NumLib
{
/// Create a IterationNumberBasedTimeStepping time stepper from the given
/// configuration.
std::unique_ptr<TimeStepAlgorithm> createIterationNumberBasedTimeStepping(
    BaseLib::ConfigTree const& config);
}  // namespace NumLib
