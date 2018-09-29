/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  \file   CreateEvolutionaryPIDcontroller.h
 *  Created on June 26, 2017, 4:43 PM
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

/// Create an EvolutionaryPIDcontroller time stepper from the given
/// configuration
std::unique_ptr<TimeStepAlgorithm> createEvolutionaryPIDcontroller(
    BaseLib::ConfigTree const& config);
}  // end of namespace NumLib
