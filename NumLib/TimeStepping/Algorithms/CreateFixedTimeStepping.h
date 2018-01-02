/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  \file   CreateFixedTimeStepping.h
 *  Created on June 26, 2017, 5:02 PM
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

/// Create a FixedTimeStepping time stepper from the given
/// configuration
std::unique_ptr<TimeStepAlgorithm> createFixedTimeStepping(
    BaseLib::ConfigTree const& config);
}  // end of namespace NumLib
