/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef NUMLIB_TIMEDISCRETIZATION_BUILDER_H
#define NUMLIB_TIMEDISCRETIZATION_BUILDER_H

#include <memory>

#include "BaseLib/ConfigTree.h"
#include "TimeDiscretization.h"

namespace NumLib
{

std::unique_ptr<TimeDiscretization> createTimeDiscretization(
    BaseLib::ConfigTree const& config);
}

#endif  // NUMLIB_TIMEDISCRETIZATION_BUILDER_H
