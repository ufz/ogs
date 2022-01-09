/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <memory>

namespace BaseLib
{
class ConfigTree;
}
namespace NumLib
{
class TimeDiscretization;
}

namespace NumLib
{
std::unique_ptr<TimeDiscretization> createTimeDiscretization(
    BaseLib::ConfigTree const& config);
}
