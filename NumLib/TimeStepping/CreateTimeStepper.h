/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  \file   CreateTimeStepper.h
 *  Created on May 2, 2017, 12:18 PM
 */

#pragma once

#include <memory>

namespace BaseLib
{
class ConfigTree;
}

namespace NumLib
{
class ITimeStepAlgorithm;
std::unique_ptr<ITimeStepAlgorithm> createTimeStepper(
    BaseLib::ConfigTree const& config);
};
